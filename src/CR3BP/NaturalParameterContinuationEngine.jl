"""
CR3BP natural parameter continuation engine wrapper

Author: Jonathan Richmond
C: 1/4/23
U: 1/27/25
"""

import MBD: CR3BPNaturalParameterContinuationEngine

export addEndCheck!, addJumpCheck!, computeFullStep, constrainNextGuess!, convergeInitialSolution
export doContinuation!, endContinuation, resetEngine!, tryConverging!

"""
    addEndCheck!(naturalParameterContinuationEngine, endCheck)

Return natural parameter continuation engine with updated end checks

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `endCheck::AbstractContinuationEndCheck`: Continuation end check
"""
function addEndCheck!(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, endCheck::MBD.AbstractContinuationEndCheck)
    push!(naturalParameterContinuationEngine.endChecks, endCheck)
end

"""
    addJumpCheck!(naturalParameterContinuationEngine, jumpCheck)

Return natural parameter continuation engine with updated jump checks

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `jumpCheck::AbstractContinuationJumpCheck`: Continuation jump check
"""
function addJumpCheck!(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, jumpCheck::MBD.AbstractContinuationJumpCheck)
    push!(naturalParameterContinuationEngine.jumpChecks, jumpCheck)
end

"""
    computeFullStep(naturalParameterContinuationEngine, data)

Return update step for free variable vector

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function computeFullStep(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, data::MBD.CR3BPContinuationData)
    fullStep::Vector{Float64} = zeros(Float64, getNumFreeVariables!(data.previousSolution))
    for (index1::MBD.Variable, value1::Int16) in data.twoPreviousSolution.freeVariableIndexMap
        for (index2::MBD.Variable, value2::Int16) in data.previousSolution.freeVariableIndexMap
            if index2.name == index1.name
                data1::Vector{Float64} = getFreeVariableData(index1)
                data2::Vector{Float64} = getFreeVariableData(index2)
                fullStep[value2:value2+length(data2)-1] = (data2-data1)./data.currentStepSize
                (index2.name == naturalParameterContinuationEngine.stepSizeGenerator.elementName) && (fullStep[value2+naturalParameterContinuationEngine.stepSizeGenerator.elementIndex-1] = 1)
            end
        end
    end

    return fullStep
end

"""
    constrainNextGuess!(naturalParameterContinuationEngine, data)

Return natural parameter continuation engine object with updated constraints

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function constrainNextGuess!(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, data::MBD.CR3BPContinuationData)
    for constraint::MBD.AbstractConstraint in keys(data.nextGuess.constraintIndexMap)
        if typeof(constraint) == MBD.CR3BPStateConstraint
            if constraint.variable.name == naturalParameterContinuationEngine.stepSizeGenerator.elementName
                variableMask::Vector{Bool} = constraint.variable.freeVariableMask
                freeCounter::Int16 = 0
                continuationIndex::Int16 = 0
                for i::Int16 in Int16(1):Int16(length(variableMask))
                    (variableMask[i] == true) && (freeCounter += 1)
                    if freeCounter == naturalParameterContinuationEngine.stepSizeGenerator.elementIndex
                        continuationIndex = i
                        break
                    end
                end
                for i::Int16 in Int16(1):Int16(length(constraint.constrainedIndices))
                    (constraint.constrainedIndices[i] == continuationIndex) && (constraint.values[i] += data.currentStepSize)
                end
            end
        end
    end
end

"""
    convergeInitialSolution(naturalParameterContinuationEngine, initialGuess)

Return converged initial solution

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `initialGuess::CR3BPMultipleShooterProblem`: Initial family member
"""
function convergeInitialSolution(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, initialGuess::MBD.CR3BPMultipleShooterProblem)
    return solve!(naturalParameterContinuationEngine.corrector, initialGuess)
end

"""
    doContinuation!(naturalParameterContinuationEngine, initialGuess1, initialGuess2)

Return family of solutions

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `initialGuess1::CR3BPMultipleShooterProblem`: First member of family
- `initialGuess2::CR3BPMultipleShooterProblem`: Second member of family
"""
function doContinuation!(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, initialGuess1::MBD.CR3BPMultipleShooterProblem, initialGuess2::MBD.CR3BPMultipleShooterProblem)
    isempty(naturalParameterContinuationEngine.endChecks) && throw(ErrorException("Cannot do continuation without at least one end check"))
    resetEngine!(naturalParameterContinuationEngine, initialGuess1, initialGuess2)
    naturalParameterContinuationEngine.corrector.printProgress = naturalParameterContinuationEngine.printProgress
    naturalParameterContinuationEngine.printProgress && println("Converging initial guesses...")
    naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution = convergeInitialSolution(naturalParameterContinuationEngine, initialGuess1)
    naturalParameterContinuationEngine.dataInProgress.previousSolution = convergeInitialSolution(naturalParameterContinuationEngine, initialGuess2)
    naturalParameterContinuationEngine.dataInProgress.numIterations = naturalParameterContinuationEngine.corrector.recentIterationCount
    push!(naturalParameterContinuationEngine.dataInProgress.family.nodes, [shallowClone(naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution.nodes[n]) for n = 1:length(naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution.nodes)], [shallowClone(naturalParameterContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(naturalParameterContinuationEngine.dataInProgress.previousSolution.nodes)])
    push!(naturalParameterContinuationEngine.dataInProgress.family.segments, [shallowClone(naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution.segments[s]) for s = 1:length(naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution.segments)], [shallowClone(naturalParameterContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(naturalParameterContinuationEngine.dataInProgress.previousSolution.segments)])
    naturalParameterContinuationEngine.dataInProgress.initialGuess = initialGuess2
    naturalParameterContinuationEngine.dataInProgress.converging = true
    naturalParameterContinuationEngine.dataInProgress.forceEndContinuation = false
    naturalParameterContinuationEngine.dataInProgress.currentStepSize = naturalParameterContinuationEngine.stepSizeGenerator.initialStepSize
    while (!endContinuation(naturalParameterContinuationEngine, naturalParameterContinuationEngine.dataInProgress) && !naturalParameterContinuationEngine.dataInProgress.forceEndContinuation)
        naturalParameterContinuationEngine.printProgress && println("\nConverging family member $(getNumSteps(naturalParameterContinuationEngine.dataInProgress)+1)...")
        naturalParameterContinuationEngine.dataInProgress.fullStep = computeFullStep(naturalParameterContinuationEngine, naturalParameterContinuationEngine.dataInProgress)
        tryConverging!(naturalParameterContinuationEngine)
        while (!naturalParameterContinuationEngine.dataInProgress.converging && !naturalParameterContinuationEngine.dataInProgress.forceEndContinuation)
            tryConverging!(naturalParameterContinuationEngine)
        end
        if (naturalParameterContinuationEngine.storeIntermediateMembers && naturalParameterContinuationEngine.dataInProgress.converging)
            push!(naturalParameterContinuationEngine.dataInProgress.family.nodes, [shallowClone(naturalParameterContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(naturalParameterContinuationEngine.dataInProgress.previousSolution.nodes)])
            push!(naturalParameterContinuationEngine.dataInProgress.family.segments, [shallowClone(naturalParameterContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(naturalParameterContinuationEngine.dataInProgress.previousSolution.segments)])
        end
    end
    if (!naturalParameterContinuationEngine.dataInProgress.converging && (getNumSteps(naturalParameterContinuationEngine.dataInProgress) == 2))
        throw(ErrorException("Could not converge any solutions beyond initial guess"))
    end
    if (!naturalParameterContinuationEngine.storeIntermediateMembers && (getNumSteps(naturalParameterContinuationEngine.dataInProgress) > 2))
        push!(naturalParameterContinuationEngine.dataInProgress.family.nodes, [shallowClone(naturalParameterContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(naturalParameterContinuationEngine.dataInProgress.previousSolution.nodes)])
            push!(naturalParameterContinuationEngine.dataInProgress.family.segments, [shallowClone(naturalParameterContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(naturalParameterContinuationEngine.dataInProgress.previousSolution.segments)])
    end

    return naturalParameterContinuationEngine.dataInProgress.family
end

"""
    endContinuation(naturalParameterContinuationEngine, data)

Return true if continuation should end

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function endContinuation(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, data::MBD.CR3BPContinuationData)
    for endCheck::MBD.AbstractContinuationEndCheck in naturalParameterContinuationEngine.endChecks
        isContinuationDone(endCheck, data) && (return true)
    end

    return false
end

"""
    resetEngine!(naturalParameterContinuationEngine, solution1, solution2)

Return natural parameter continuation engine with reset data

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
- `solution1::CR3BPMultipleShooterProblem`: First member of family
- `solution2::CR3BPMultipleShooterProblem`: Second member of family
"""
function resetEngine!(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine, solution1::MBD.CR3BPMultipleShooterProblem, solution2::MBD.CR3BPMultipleShooterProblem)
    naturalParameterContinuationEngine.dataInProgress = MBD.CR3BPContinuationData(solution1, solution2)
end

"""
    tryConverging!(naturalParameterContinuationEngine, solution)

Return updated natural parameter continuation engine object

# Arguments
- `naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine`: CR3BP natural parameter continuation engine object
"""
function tryConverging!(naturalParameterContinuationEngine::CR3BPNaturalParameterContinuationEngine)
    updateStepSize!(naturalParameterContinuationEngine.stepSizeGenerator, naturalParameterContinuationEngine.dataInProgress)
    naturalParameterContinuationEngine.printProgress && println("Current step size: $(naturalParameterContinuationEngine.dataInProgress.currentStepSize)")
    naturalParameterContinuationEngine.dataInProgress.nextGuess = deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution)
    setFreeVariableVector!(naturalParameterContinuationEngine.dataInProgress.nextGuess, getFreeVariableVector!(naturalParameterContinuationEngine.dataInProgress.previousSolution)+naturalParameterContinuationEngine.dataInProgress.fullStep.*naturalParameterContinuationEngine.dataInProgress.currentStepSize)
    constrainNextGuess!(naturalParameterContinuationEngine, naturalParameterContinuationEngine.dataInProgress)
    try
        twoPreviousConvergedSolution::MBD.CR3BPMultipleShooterProblem = naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution
        previousConvergedSolution::MBD.CR3BPMultipleShooterProblem = naturalParameterContinuationEngine.dataInProgress.previousSolution
        naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution = deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution)
        naturalParameterContinuationEngine.dataInProgress.previousSolution = solve!(naturalParameterContinuationEngine.corrector, naturalParameterContinuationEngine.dataInProgress.nextGuess)
        naturalParameterContinuationEngine.dataInProgress.converging = true
        for jumpCheck::MBD.AbstractContinuationJumpCheck in naturalParameterContinuationEngine.jumpChecks
            if typeof(jumpCheck) == MBD.BoundingBoxJumpCheck
                for (index::MBD.Variable, value::Int16) in naturalParameterContinuationEngine.dataInProgress.previousSolution.freeVariableIndexMap
                    if index.name == jumpCheck.paramName
                        addBounds!(jumpCheck, naturalParameterContinuationEngine.dataInProgress.previousSolution, index, jumpCheck.paramBounds)
                        naturalParameterContinuationEngine.dataInProgress.converging = isFamilyMember(jumpCheck, naturalParameterContinuationEngine.dataInProgress)
                        !naturalParameterContinuationEngine.dataInProgress.converging && println("Solution jumped")
                        removeBounds!(jumpCheck, naturalParameterContinuationEngine.dataInProgress.previousSolution, index)
                    end
                end
            end
        end
        if naturalParameterContinuationEngine.dataInProgress.converging
            naturalParameterContinuationEngine.dataInProgress.numIterations = naturalParameterContinuationEngine.corrector.recentIterationCount
        else
            naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution = twoPreviousConvergedSolution
            naturalParameterContinuationEngine.dataInProgress.previousSolution = previousConvergedSolution
        end
    catch
        naturalParameterContinuationEngine.dataInProgress.converging = false
        println("Failed to converge")
    end
end
