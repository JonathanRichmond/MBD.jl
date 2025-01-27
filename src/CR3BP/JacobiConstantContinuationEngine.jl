"""
Jacobi constant continuation engine wrapper

Author: Jonathan Richmond
C: 1/11/23
U: 1/26/25
"""

import MBD: JacobiConstantContinuationEngine

export addEndCheck!, addJumpCheck!, computeFullStep, constrainNextGuess!, convergeInitialSolution
export doContinuation!, endContinuation, resetEngine!, tryConverging!

"""
    addEndCheck!(jacobiConstantContinuationEngine, endCheck)

Return Jacobi constant continuation engine with updated end checks

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `endCheck::AbstractContinuationEndCheck`: Continuation end check
"""
function addEndCheck!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, endCheck::MBD.AbstractContinuationEndCheck)
    push!(jacobiConstantContinuationEngine.endChecks, endCheck)
end

"""
    addJumpCheck!(jacobiConstantContinuationEngine, jumpCheck)

Return Jacobi constant continuation engine with updated jump checks

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `jumpCheck::AbstractContinuationJumpCheck`: Continuation jump check
"""
function addJumpCheck!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, jumpCheck::MBD.AbstractContinuationJumpCheck)
    push!(jacobiConstantContinuationEngine.jumpChecks, jumpCheck)
end

"""
    computeFullStep(jacobiConstantContinuationEngine, data)

Return update step for free variable vector

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function computeFullStep(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, data::MBD.CR3BPContinuationData)
    fullStep::Vector{Float64} = zeros(Float64, getNumFreeVariables!(data.previousSolution))
    for (index1::MBD.Variable, value1::Int16) in data.twoPreviousSolution.freeVariableIndexMap
        for (index2::MBD.Variable, value2::Int16) in data.previousSolution.freeVariableIndexMap
            if index2.name == index1.name
                data1::Vector{Float64} = getFreeVariableData(index1)
                data2::Vector{Float64} = getFreeVariableData(index2)
                fullStep[value2:value2+length(data2)-1] = (data2-data1)./data.currentStepSize
            end
        end
    end

    return fullStep
end

"""
    constrainNextGuess!(jacobiConstantContinuationEngine, data)

Return Jacobi constant continuation engine object with updated constraints

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function constrainNextGuess!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, data::MBD.CR3BPContinuationData)
    for constraint::MBD.AbstractConstraint in keys(data.nextGuess.constraintIndexMap)
        if typeof(constraint) == MBD.JacobiConstraint
            constraint.value += data.currentStepSize
        end
    end
end

"""
    convergeInitialSolution(jacobiConstantContinuationEngine, initialGuess)

Return converged initial solution

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `initialGuess::CR3BPMultipleShooterProblem`: Initial family member
"""
function convergeInitialSolution(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, initialGuess::MBD.CR3BPMultipleShooterProblem)
    return solve!(jacobiConstantContinuationEngine.corrector, initialGuess)
end

"""
    doContinuation!(jacobiConstantContinuationEngine, initialGuess1, initialGuess2)

Return family of solutions

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `initialGuess1::CR3BPMultipleShooterProblem`: First member of family
- `initialGuess2::CR3BPMultipleShooterProblem`: Second member of family
"""
function doContinuation!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, initialGuess1::MBD.CR3BPMultipleShooterProblem, initialGuess2::MBD.CR3BPMultipleShooterProblem)
    isempty(jacobiConstantContinuationEngine.endChecks) && throw(ErrorException("Cannot do continuation without at least one end check"))
    resetEngine!(jacobiConstantContinuationEngine, initialGuess1, initialGuess2)
    jacobiConstantContinuationEngine.corrector.printProgress = jacobiConstantContinuationEngine.printProgress
    jacobiConstantContinuationEngine.printProgress && println("Converging initial guesses...")
    jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution = convergeInitialSolution(jacobiConstantContinuationEngine, initialGuess1)
    jacobiConstantContinuationEngine.dataInProgress.previousSolution = convergeInitialSolution(jacobiConstantContinuationEngine, initialGuess2)
    jacobiConstantContinuationEngine.dataInProgress.numIterations = jacobiConstantContinuationEngine.corrector.recentIterationCount
    push!(jacobiConstantContinuationEngine.dataInProgress.family.nodes, [shallowClone(jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution.nodes[n]) for n = 1:length(jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution.nodes)], [shallowClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(jacobiConstantContinuationEngine.dataInProgress.previousSolution.nodes)])
    push!(jacobiConstantContinuationEngine.dataInProgress.family.segments, [shallowClone(jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution.segments[s]) for s = 1:length(jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution.segments)], [shallowClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(jacobiConstantContinuationEngine.dataInProgress.previousSolution.segments)])
    jacobiConstantContinuationEngine.dataInProgress.initialGuess = initialGuess2
    jacobiConstantContinuationEngine.dataInProgress.converging = true
    jacobiConstantContinuationEngine.dataInProgress.forceEndContinuation = false
    jacobiConstantContinuationEngine.dataInProgress.currentStepSize = jacobiConstantContinuationEngine.stepSizeGenerator.initialStepSize
    while (!endContinuation(jacobiConstantContinuationEngine, jacobiConstantContinuationEngine.dataInProgress) && !jacobiConstantContinuationEngine.dataInProgress.forceEndContinuation)
        jacobiConstantContinuationEngine.printProgress && println("\nConverging family member $(getNumSteps(jacobiConstantContinuationEngine.dataInProgress)+1)...")
        jacobiConstantContinuationEngine.dataInProgress.fullStep = computeFullStep(jacobiConstantContinuationEngine, jacobiConstantContinuationEngine.dataInProgress)
        tryConverging!(jacobiConstantContinuationEngine)
        while (!jacobiConstantContinuationEngine.dataInProgress.converging && !jacobiConstantContinuationEngine.dataInProgress.forceEndContinuation)
            tryConverging!(jacobiConstantContinuationEngine)
        end
        if (jacobiConstantContinuationEngine.storeIntermediateMembers && jacobiConstantContinuationEngine.dataInProgress.converging)
            push!(jacobiConstantContinuationEngine.dataInProgress.family.nodes, [shallowClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(jacobiConstantContinuationEngine.dataInProgress.previousSolution.nodes)])
            push!(jacobiConstantContinuationEngine.dataInProgress.family.segments, [shallowClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(jacobiConstantContinuationEngine.dataInProgress.previousSolution.segments)])
        end
    end
    if (!jacobiConstantContinuationEngine.dataInProgress.converging && (getNumSteps(jacobiConstantContinuationEngine.dataInProgress) == 2))
        throw(ErrorException("Could not converge any solutions beyond initial guess"))
    end
    if (!jacobiConstantContinuationEngine.storeIntermediateMembers && (getNumSteps(jacobiConstantContinuationEngine.dataInProgress) > 2))
        push!(jacobiConstantContinuationEngine.dataInProgress.family.nodes, [shallowClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(jacobiConstantContinuationEngine.dataInProgress.previousSolution.nodes)])
            push!(jacobiConstantContinuationEngine.dataInProgress.family.segments, [shallowClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(jacobiConstantContinuationEngine.dataInProgress.previousSolution.segments)])
    end

    return jacobiConstantContinuationEngine.dataInProgress.family
end

"""
    endContinuation(jacobiConstantContinuationEngine, data)

Return true if continuation should end

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function endContinuation(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, data::MBD.CR3BPContinuationData)
    for endCheck::MBD.AbstractContinuationEndCheck in jacobiConstantContinuationEngine.endChecks
        isContinuationDone(endCheck, data) && (return true)
    end

    return false
end

"""
    resetEngine!(jacobiConstantContinuationEngine, solution1, solution2)

Return Jacobi constant continuation engine with reset data

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `solution1::CR3BPMultipleShooterProblem`: First member of family
- `solution2::CR3BPMultipleShooterProblem`: Second member of family
"""
function resetEngine!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, solution1::MBD.CR3BPMultipleShooterProblem, solution2::MBD.CR3BPMultipleShooterProblem)
    jacobiConstantContinuationEngine.dataInProgress = MBD.CR3BPContinuationData(solution1, solution2)
end

"""
    tryConverging!(jacobiConstantContinuationEngine)

Return updated Jacobi constant continuation engine object

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
"""
function tryConverging!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine)
    updateStepSize!(jacobiConstantContinuationEngine.stepSizeGenerator, jacobiConstantContinuationEngine.dataInProgress)
    jacobiConstantContinuationEngine.printProgress && println("Current step size: $(jacobiConstantContinuationEngine.dataInProgress.currentStepSize)")
    jacobiConstantContinuationEngine.dataInProgress.nextGuess = deepClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution)
    setFreeVariableVector!(jacobiConstantContinuationEngine.dataInProgress.nextGuess, getFreeVariableVector!(jacobiConstantContinuationEngine.dataInProgress.previousSolution)+jacobiConstantContinuationEngine.dataInProgress.fullStep.*jacobiConstantContinuationEngine.dataInProgress.currentStepSize)
    constrainNextGuess!(jacobiConstantContinuationEngine, jacobiConstantContinuationEngine.dataInProgress)
    try
        twoPreviousConvergedSolution::MBD.CR3BPMultipleShooterProblem = jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution
        previousConvergedSolution::MBD.CR3BPMultipleShooterProblem = jacobiConstantContinuationEngine.dataInProgress.previousSolution
        jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution = deepClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution)
        jacobiConstantContinuationEngine.dataInProgress.previousSolution = solve!(jacobiConstantContinuationEngine.corrector, jacobiConstantContinuationEngine.dataInProgress.nextGuess)
        jacobiConstantContinuationEngine.dataInProgress.converging = true
        for jumpCheck::MBD.AbstractContinuationJumpCheck in jacobiConstantContinuationEngine.jumpChecks
            if typeof(jumpCheck) == MBD.BoundingBoxJumpCheck
                for (index::MBD.Variable, value::Int16) in jacobiConstantContinuationEngine.dataInProgress.previousSolution.freeVariableIndexMap
                    if index.name == jumpCheck.paramName
                        addBounds!(jumpCheck, jacobiConstantContinuationEngine.dataInProgress.previousSolution, index, jumpCheck.paramBounds)
                        jacobiConstantContinuationEngine.dataInProgress.converging = isFamilyMember(jumpCheck, jacobiConstantContinuationEngine.dataInProgress)
                        !jacobiConstantContinuationEngine.dataInProgress.converging && println("Solution jumped")
                        removeBounds!(jumpCheck, jacobiConstantContinuationEngine.dataInProgress.previousSolution, index)
                    end
                end
            end
        end
        if jacobiConstantContinuationEngine.dataInProgress.converging
            jacobiConstantContinuationEngine.dataInProgress.numIterations = jacobiConstantContinuationEngine.corrector.recentIterationCount
        else
            jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution = twoPreviousConvergedSolution
            jacobiConstantContinuationEngine.dataInProgress.previousSolution = previousConvergedSolution
        end
    catch
        jacobiConstantContinuationEngine.dataInProgress.converging = false
        println("Failed to converge")
    end
end
