"""
Natural parameter continuation engine wrapper

Author: Jonathan Richmond
C: 1/4/23
U: 9/8/23
"""

import MBD: NaturalParameterContinuationEngine

export addEndCheck!, addJumpCheck!, computeFullStep, constrainNextGuess!
export convergeInitialSolution, doContinuation!, endContinuation, resetEngine!
export tryConverging!

"""
    addEndCheck!(naturalParameterContinuationEngine, endCheck)

Return natural parameter continuation engine with updated end checks

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `endCheck::AbstractContinuationEndCheck`: Continuation end check
"""
function addEndCheck!(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, endCheck::MBD.AbstractContinuationEndCheck)
    push!(naturalParameterContinuationEngine.endChecks, endCheck)
end

"""
    addJumpCheck!(naturalParameterContinuationEngine, jumpCheck)

Return natural parameter continuation engine with updated jump checks

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `jumpCheck::AbstractContinuationJumpCheck`: Continuation jump check
"""
function addJumpCheck!(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, jumpCheck::MBD.AbstractContinuationJumpCheck)
    push!(naturalParameterContinuationEngine.jumpChecks, jumpCheck)
end

"""
    computeFullStep(naturalParameterContinuationEngine, data)

Return update step for free variable vector

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `data::ContinuationData`: Continuation data object
"""
function computeFullStep(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, data::MBD.ContinuationData)
    fullStep::Vector{Float64} = zeros(Float64, getNumberFreeVariables(data.previousSolution))
    for (index1::MBD.Variable, value1::Int64) in data.twoPreviousSolution.freeVariableIndexMap
        for (index2::MBD.Variable, value2::Int64) in data.previousSolution.freeVariableIndexMap
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
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `data::ContinuationData`: Continuation data object
"""
function constrainNextGuess!(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, data::MBD.ContinuationData)
    for constraint::MBD.AbstractConstraint in keys(data.nextGuess.constraintIndexMap)
        if typeof(constraint) == MBD.StateConstraint
            if constraint.variable.name == naturalParameterContinuationEngine.stepSizeGenerator.elementName
                variableMask::Vector{Bool} = constraint.variable.freeVariableMask
                freeCounter::Int64 = 0
                continuationIndex::Int64 = 0
                for i::Int64 in 1:length(variableMask)
                    (variableMask[i] == true) && (freeCounter += 1)
                    (freeCounter == naturalParameterContinuationEngine.stepSizeGenerator.elementIndex) && (continuationIndex = i)
                end
                for i::Int64 in 1:length(constraint.constrainedIndices)
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
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `initialGuess::MultipleShooterProblem`: Initial family member
"""
function convergeInitialSolution(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, initialGuess::MBD.MultipleShooterProblem)
    return solve!(naturalParameterContinuationEngine.corrector, initialGuess)
end

"""
    doContinuation!(naturalParameterContinuationEngine, initialGuess1, initialGuess2)

Return family of solutions

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `initialGuess1::MultipleShooterProblem`: First member of family
- `initialGuess2::MultipleShooterProblem`: Second member of family
"""
function doContinuation!(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, initialGuess1::MBD.MultipleShooterProblem, initialGuess2::MBD.MultipleShooterProblem)
    isempty(naturalParameterContinuationEngine.endChecks) && throw(ErrorException("Cannot do continuation without at least one end check"))
    resetEngine!(naturalParameterContinuationEngine)
    setPrintProgress!(naturalParameterContinuationEngine.corrector, naturalParameterContinuationEngine.printProgress)
    naturalParameterContinuationEngine.printProgress && println("Converging initial guesses...")
    naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution = convergeInitialSolution(naturalParameterContinuationEngine, initialGuess1)
    naturalParameterContinuationEngine.dataInProgress.previousSolution = convergeInitialSolution(naturalParameterContinuationEngine, initialGuess2)
    naturalParameterContinuationEngine.dataInProgress.numIterations = naturalParameterContinuationEngine.corrector.recentIterationCount
    push!(naturalParameterContinuationEngine.dataInProgress.familyMembers, deepClone(naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution), deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution))
    naturalParameterContinuationEngine.dataInProgress.initialGuess = initialGuess2
    naturalParameterContinuationEngine.dataInProgress.stepCount = 2
    naturalParameterContinuationEngine.dataInProgress.converging = true
    naturalParameterContinuationEngine.dataInProgress.forceEndContinuation = false
    naturalParameterContinuationEngine.dataInProgress.currentStepSize = naturalParameterContinuationEngine.stepSizeGenerator.initialStepSize
    while (!endContinuation(naturalParameterContinuationEngine, naturalParameterContinuationEngine.dataInProgress) && !naturalParameterContinuationEngine.dataInProgress.forceEndContinuation)
        naturalParameterContinuationEngine.printProgress && println("\nConverging family member $(naturalParameterContinuationEngine.dataInProgress.stepCount+1)...")
        naturalParameterContinuationEngine.dataInProgress.fullStep = computeFullStep(naturalParameterContinuationEngine, naturalParameterContinuationEngine.dataInProgress)
        tryConverging!(naturalParameterContinuationEngine)
        while (!naturalParameterContinuationEngine.dataInProgress.converging && !naturalParameterContinuationEngine.dataInProgress.forceEndContinuation)
            tryConverging!(naturalParameterContinuationEngine)
        end
        (naturalParameterContinuationEngine.storeIntermediateMembers && naturalParameterContinuationEngine.dataInProgress.converging) && push!(naturalParameterContinuationEngine.dataInProgress.familyMembers, deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution))
    end
    if (!naturalParameterContinuationEngine.dataInProgress.converging && (naturalParameterContinuationEngine.dataInProgress.stepCount == 2))
        throw(ErrorException("Could not converge any solutions beyond initial guess"))
    end
    if (!naturalParameterContinuationEngine.storeIntermediateMembers && (naturalParameterContinuationEngine.dataInProgress.stepCount > 2))
        push!(naturalParameterContinuationEngine.dataInProgress.familyMembers, deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution))
    end

    return naturalParameterContinuationEngine.dataInProgress.familyMembers
end

"""
    endContinuation(naturalParameterContinuationEngine, data)

Return true if continuation should end

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
- `data::ContinuationData`: Continuation data object
"""
function endContinuation(naturalParameterContinuationEngine::NaturalParameterContinuationEngine, data::MBD.ContinuationData)
    for endCheck::MBD.AbstractContinuationEndCheck in naturalParameterContinuationEngine.endChecks
        isContinuationDone(endCheck, data) && (return true)
    end

    return false
end

"""
    resetEngine!(naturalParameterContinuationEngine)

Return natural parameter continuation engine with reset data

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
"""
function resetEngine!(naturalParameterContinuationEngine::NaturalParameterContinuationEngine)
    naturalParameterContinuationEngine.dataInProgress = MBD.ContinuationData()
end

"""
    tryConverging!(naturalParameterContinuationEngine)

Return updated natural parameter continuation engine object

# Arguments
- `naturalParameterContinuationEngine::NaturalParameterContinuationEngine`: Natural parameter continuation engine object
"""
function tryConverging!(naturalParameterContinuationEngine::NaturalParameterContinuationEngine)
    updateStepSize!(naturalParameterContinuationEngine.stepSizeGenerator, naturalParameterContinuationEngine.dataInProgress)
    naturalParameterContinuationEngine.printProgress && println("Current step size: $(naturalParameterContinuationEngine.dataInProgress.currentStepSize)")
    naturalParameterContinuationEngine.dataInProgress.nextGuess = deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution)
    setFreeVariableVector!(naturalParameterContinuationEngine.dataInProgress.nextGuess, getFreeVariableVector!(naturalParameterContinuationEngine.dataInProgress.previousSolution)+naturalParameterContinuationEngine.dataInProgress.fullStep.*naturalParameterContinuationEngine.dataInProgress.currentStepSize)
    constrainNextGuess!(naturalParameterContinuationEngine, naturalParameterContinuationEngine.dataInProgress)
    try
        twoPreviousConvergedSolution::MBD.MultipleShooterProblem = naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution
        previousConvergedSolution::MBD.MultipleShooterProblem = naturalParameterContinuationEngine.dataInProgress.previousSolution
        naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution = deepClone(naturalParameterContinuationEngine.dataInProgress.previousSolution)
        naturalParameterContinuationEngine.dataInProgress.previousSolution = solve!(naturalParameterContinuationEngine.corrector, naturalParameterContinuationEngine.dataInProgress.nextGuess)
        naturalParameterContinuationEngine.dataInProgress.converging = true
        for jumpCheck::MBD.AbstractContinuationJumpCheck in naturalParameterContinuationEngine.jumpChecks
            if typeof(jumpCheck) == MBD.BoundingBoxJumpCheck
                for (index::MBD.Variable, value::Int64) in naturalParameterContinuationEngine.dataInProgress.previousSolution.freeVariableIndexMap
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
            naturalParameterContinuationEngine.dataInProgress.stepCount += 1
        else
            naturalParameterContinuationEngine.dataInProgress.twoPreviousSolution = twoPreviousConvergedSolution
            naturalParameterContinuationEngine.dataInProgress.previousSolution = previousConvergedSolution
        end
    catch
        naturalParameterContinuationEngine.dataInProgress.converging = false
        println("Failed to converge")
    end
end
