"""
Jacobi constant continuation engine wrapper

Author: Jonathan Richmond
C: 1/11/23
U: 9/10/23
"""

import MBD: JacobiConstantContinuationEngine

export addEndCheck!, addJumpCheck!, computeFullStep, constrainNextGuess!
export convergeInitialSolution, doContinuation!, endContinuation, resetEngine!
export tryConverging!

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
- `data::ContinuationData`: Continuation data object
"""
function computeFullStep(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, data::MBD.ContinuationData)
    fullStep::Vector{Float64} = zeros(Float64, getNumberFreeVariables(data.previousSolution))
    for (index1::MBD.Variable, value1::Int64) in data.twoPreviousSolution.freeVariableIndexMap
        for (index2::MBD.Variable, value2::Int64) in data.previousSolution.freeVariableIndexMap
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
- `data::ContinuationData`: Continuation data object
"""
function constrainNextGuess!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, data::MBD.ContinuationData)
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
- `initialGuess::MultipleShooterProblem`: Initial family member
"""
function convergeInitialSolution(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, initialGuess::MBD.MultipleShooterProblem)
    return solve!(jacobiConstantContinuationEngine.corrector, initialGuess)
end

"""
    doContinuation!(jacobiConstantContinuationEngine, initialGuess1, initialGuess2)

Return family of solutions

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `initialGuess1::MultipleShooterProblem`: First member of family
- `initialGuess2::MultipleShooterProblem`: Second member of family
"""
function doContinuation!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, initialGuess1::MBD.MultipleShooterProblem, initialGuess2::MBD.MultipleShooterProblem)
    isempty(jacobiConstantContinuationEngine.endChecks) && throw(ErrorException("Cannot do continuation without at least one end check"))
    resetEngine!(jacobiConstantContinuationEngine)
    setPrintProgress!(jacobiConstantContinuationEngine.corrector, jacobiConstantContinuationEngine.printProgress)
    jacobiConstantContinuationEngine.printProgress && println("Converging initial guesses...")
    jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution = convergeInitialSolution(jacobiConstantContinuationEngine, initialGuess1)
    jacobiConstantContinuationEngine.dataInProgress.previousSolution = convergeInitialSolution(jacobiConstantContinuationEngine, initialGuess2)
    jacobiConstantContinuationEngine.dataInProgress.numIterations = jacobiConstantContinuationEngine.corrector.recentIterationCount
    push!(jacobiConstantContinuationEngine.dataInProgress.familyMembers, deepClone(jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution), deepClone(jacobiConstantContinuationEngine.previousSolution))
    jacobiConstantContinuationEngine.dataInProgress.initialGuess = initialGuess2
    jacobiConstantContinuationEngine.dataInProgress.stepCount = 2
    jacobiConstantContinuationEngine.dataInProgress.converging = true
    jacobiConstantContinuationEngine.dataInProgress.forceEndContinuation = false
    jacobiConstantContinuationEngine.dataInProgress.currentStepSize = jacobiConstantContinuationEngine.stepSizeGenerator.initialStepSize
    while (!endContinuation(jacobiConstantContinuationEngine, jacobiConstantContinuationEngine.dataInProgress) && !jacobiConstantContinuationEngine.dataInProgress.forceEndContinuation)
        jacobiConstantContinuationEngine.printProgress && println("\nConverging family member $(jacobiConstantContinuationEngine.dataInProgress.stepCount+1)...")
        jacobiConstantContinuationEngine.dataInProgress.fullStep = computeFullStep(jacobiConstantContinuationEngine, jacobiConstantContinuationEngine.dataInProgress)
        tryConverging!(jacobiConstantContinuationEngine)
        while (!jacobiConstantContinuationEngine.dataInProgress.converging && !jacobiConstantContinuationEngine.dataInProgress.forceEndContinuation)
            tryConverging!(jacobiConstantContinuationEngine)
        end
        (jacobiConstantContinuationEngine.storeIntermediateMembers && jacobiConstantContinuationEngine.dataInProgress.converging) && push!(jacobiConstantContinuationEngine.dataInProgress.familyMembers, deepClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution))
    end
    if (!jacobiConstantContinuationEngine.dataInProgress.converging && (jacobiConstantContinuationEngine.dataInProgress.stepCount == 2))
        throw(ErrorException("Could not converge any solutions beyond initial guess"))
    end
    if (!jacobiConstantContinuationEngine.storeIntermediateMembers && (jacobiConstantContinuationEngine.dataInProgress.stepCount > 2))
        push!(jacobiConstantContinuationEngine.dataInProgress.familyMembers, deepClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution))
    end

    return jacobiConstantContinuationEngine.dataInProgress.familyMembers
end

"""
    endContinuation(jacobiConstantContinuationEngine, data)

Return true if continuation should end

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
- `data::ContinuationData`: Continuation data object
"""
function endContinuation(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine, data::MBD.ContinuationData)
    for endCheck::MBD.AbstractContinuationEndCheck in jacobiConstantContinuationEngine.endChecks
        isContinuationDone(endCheck, data) && (return true)
    end

    return false
end

"""
    resetEngine!(jacobiConstantContinuationEngine)

Return Jacobi constant continuation engine with reset data

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
"""
function resetEngine!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine)
    jacobiConstantContinuationEngine.dataInProgress = ContinuationData()
end

"""
    tryConverging!(jacobiConstantContinuationEngine)

Return updated Jacobi constant continuation engine object

# Arguments
- `jacobiConstantContinuationEngine::JacobiConstantContinuationEngine`: Jacobi constant continuation engine object
"""
function tryConverging!(jacobiConstantContinuationEngine::JacobiConstantContinuationEngine)
    updateStepSize!(jacobiConstantContinuationEngine.stepSizeGenerator, jacobiConstantContinuationEngine.dataInProgress)
    println("Current step size: $(jacobiConstantContinuationEngine.dataInProgress.currentStepSize)")
    jacobiConstantContinuationEngine.dataInProgress.nextGuess = deepClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution)
    setFreeVariableVector!(jacobiConstantContinuationEngine.dataInProgress.nextGuess, getFreeVariableVector!(jacobiConstantContinuationEngine.dataInProgress.previousSolution)+jacobiConstantContinuationEngine.dataInProgress.fullStep.*jacobiConstantContinuationEngine.dataInProgress.currentStepSize)
    constrainNextGuess!(jacobiConstantContinuationEngine, jacobiConstantContinuationEngine.dataInProgress)
    try
        twoPreviousConvergedSolution::MBD.MultipleShooterProblem = jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution
        previousConvergedSolution::MBD.MultipleShooterProblem = jacobiConstantContinuationEngine.dataInProgress.previousSolution
        jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution = deepClone(jacobiConstantContinuationEngine.dataInProgress.previousSolution)
        jacobiConstantContinuationEngine.dataInProgress.previousSolution = solve!(jacobiConstantContinuationEngine.corrector, jacobiConstantContinuationEngine.dataInProgress.nextGuess)
        jacobiConstantContinuationEngine.dataInProgress.converging = true
        for jumpCheck::MBD.AbstractContinuationJumpCheck in jacobiConstantContinuationEngine.jumpChecks
            if typeof(jumpCheck) == MBD.BoundingBoxJumpCheck
                for (index::MBD.Variable, value::Int64) in jacobiConstantContinuationEngine.dataInProgress.previousSolution.freeVariableIndexMap
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
            jacobiConstantContinuationEngine.dataInProgress.stepCount += 1
        else
            jacobiConstantContinuationEngine.dataInProgress.twoPreviousSolution = twoPreviousConvergedSolution
            jacobiConstantContinuationEngine.dataInProgress.previousSolution = previousConvergedSolution
        end
    catch
        jacobiConstantContinuationEngine.dataInProgress.converging = false
        println("Failed to converge")
    end
end
