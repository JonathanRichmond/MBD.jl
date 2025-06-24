"""
P4 mass continuation engine wrapper

Author: Jonathan Richmond
C: 6/18/23
"""

import LinearAlgebra
import MBD: P4MassContinuationEngine

export addEndCheck!, addJumpCheck!, computeFullStep, constrainNextGuess!, convergeInitialSolution
export doContinuation!, endContinuation, resetEngine!, tryConverging!

"""
    addEndCheck!(p4MassContinuationEngine, endCheck)

Return P4 mass continuation engine with updated end checks

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `endCheck::AbstractContinuationEndCheck`: Continuation end check
"""
function addEndCheck!(p4MassContinuationEngine::P4MassContinuationEngine, endCheck::MBD.AbstractContinuationEndCheck)
    push!(p4MassContinuationEngine.endChecks, endCheck)
end

"""
    addJumpCheck!(p4MassContinuationEngine, jumpCheck)

Return P4 mass continuation engine with updated jump checks

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `jumpCheck::AbstractContinuationJumpCheck`: Continuation jump check
"""
function addJumpCheck!(p4MassContinuationEngine::P4MassContinuationEngine, jumpCheck::MBD.AbstractContinuationJumpCheck)
    push!(p4MassContinuationEngine.jumpChecks, jumpCheck)
    println("Jump check added")
end

"""
    computeFullStep(p4MassContinuationEngine, data)

Return update step for free variable vector

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `data::BCR4BP12ContinuationData`: BCR4BP P1-P2 continuation data object
"""
function computeFullStep(p4MassContinuationEngine::P4MassContinuationEngine, data::MBD.BCR4BP12ContinuationData)
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
    constrainNextGuess!(p4MassContinuationEngine, data)

Return P4 mass continuation engine object with updated constraints

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `data::BCR4BP12ContinuationData`: BCR4BP P1-P2 continuation data object
"""
function constrainNextGuess!(p4MassContinuationEngine::P4MassContinuationEngine, data::MBD.BCR4BP12ContinuationData)
    for node::MBD.BCR4BP12Node in data.nextGuess.nodes
        node.dynamicsModel.systemData.P4Mass += data.currentStepSize*node.dynamicsModel.systemData.primaryData[3].mass
    end
end

"""
    convergeInitialSolution(p4MassContinuationEngine, initialGuess)

Return converged initial solution

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `initialGuess::BCR4BP12MultipleShooterProblem`: Initial family member
"""
function convergeInitialSolution(p4MassContinuationEngine::P4MassContinuationEngine, initialGuess::MBD.BCR4BP12MultipleShooterProblem)
    return solve!(p4MassContinuationEngine.corrector, initialGuess)
end

"""
    doContinuation!(p4MassContinuationEngine, initialGuess1, initialGuess2)

Return family of solutions

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `initialGuess1::BCR4BP12MultipleShooterProblem`: First member of family
- `initialGuess2::BCR4BP12MultipleShooterProblem`: Second member of family
"""
function doContinuation!(p4MassContinuationEngine::P4MassContinuationEngine, initialGuess1::MBD.BCR4BP12MultipleShooterProblem, initialGuess2::MBD.BCR4BP12MultipleShooterProblem)
    isempty(p4MassContinuationEngine.endChecks) && throw(ErrorException("Cannot do continuation without at least one end check"))
    resetEngine!(p4MassContinuationEngine, initialGuess1, initialGuess2)
    p4MassContinuationEngine.corrector.printProgress = p4MassContinuationEngine.printProgress
    p4MassContinuationEngine.printProgress && println("Converging initial guesses...")
    p4MassContinuationEngine.dataInProgress.twoPreviousSolution = convergeInitialSolution(p4MassContinuationEngine, initialGuess1)
    p4MassContinuationEngine.dataInProgress.previousSolution = convergeInitialSolution(p4MassContinuationEngine, initialGuess2)
    p4MassContinuationEngine.dataInProgress.numIterations = p4MassContinuationEngine.corrector.recentIterationCount
    push!(p4MassContinuationEngine.dataInProgress.family.nodes, [shallowClone(p4MassContinuationEngine.dataInProgress.twoPreviousSolution.nodes[n]) for n = 1:length(p4MassContinuationEngine.dataInProgress.twoPreviousSolution.nodes)], [shallowClone(p4MassContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(p4MassContinuationEngine.dataInProgress.previousSolution.nodes)])
    push!(p4MassContinuationEngine.dataInProgress.family.segments, [shallowClone(p4MassContinuationEngine.dataInProgress.twoPreviousSolution.segments[s]) for s = 1:length(p4MassContinuationEngine.dataInProgress.twoPreviousSolution.segments)], [shallowClone(p4MassContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(p4MassContinuationEngine.dataInProgress.previousSolution.segments)])
    p4MassContinuationEngine.dataInProgress.initialGuess = initialGuess2
    p4MassContinuationEngine.dataInProgress.converging = true
    p4MassContinuationEngine.dataInProgress.forceEndContinuation = false
    p4MassContinuationEngine.dataInProgress.currentStepSize = p4MassContinuationEngine.stepSizeGenerator.initialStepSize
    while (!endContinuation(p4MassContinuationEngine, p4MassContinuationEngine.dataInProgress) && !p4MassContinuationEngine.dataInProgress.forceEndContinuation)
        p4MassContinuationEngine.printProgress && println("\nConverging family member $(getNumSteps(p4MassContinuationEngine.dataInProgress)+1)...")
        p4MassContinuationEngine.dataInProgress.fullStep = computeFullStep(p4MassContinuationEngine, p4MassContinuationEngine.dataInProgress)
        tryConverging!(p4MassContinuationEngine)
        while (!p4MassContinuationEngine.dataInProgress.converging && !p4MassContinuationEngine.dataInProgress.forceEndContinuation)
            tryConverging!(p4MassContinuationEngine)
        end
        if (p4MassContinuationEngine.storeIntermediateMembers && p4MassContinuationEngine.dataInProgress.converging)
            push!(p4MassContinuationEngine.dataInProgress.family.nodes, [shallowClone(p4MassContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(p4MassContinuationEngine.dataInProgress.previousSolution.nodes)])
            push!(p4MassContinuationEngine.dataInProgress.family.segments, [shallowClone(p4MassContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(p4MassContinuationEngine.dataInProgress.previousSolution.segments)])
        end
    end
    if (!p4MassContinuationEngine.dataInProgress.converging && (getNumSteps(p4MassContinuationEngine.dataInProgress) == 2))
        throw(ErrorException("Could not converge any solutions beyond initial guess"))
    end
    if (!p4MassContinuationEngine.storeIntermediateMembers && (getNumSteps(p4MassContinuationEngine.dataInProgress) > 2))
        push!(p4MassContinuationEngine.dataInProgress.family.nodes, [shallowClone(p4MassContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(p4MassContinuationEngine.dataInProgress.previousSolution.nodes)])
        push!(p4MassContinuationEngine.dataInProgress.family.segments, [shallowClone(p4MassContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(p4MassContinuationEngine.dataInProgress.previousSolution.segments)])
    end

    return p4MassContinuationEngine.dataInProgress.family
end

"""
    endContinuation(p4MassContinuationEngine, data)

Return true if continuation should end

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `data::BCR4BP12ContinuationData`: BCR4BP P1-P2 continuation data object
"""
function endContinuation(p4MassContinuationEngine::P4MassContinuationEngine, data::MBD.BCR4BP12ContinuationData)
    for endCheck::MBD.AbstractContinuationEndCheck in p4MassContinuationEngine.endChecks
        isContinuationDone(endCheck, data) && (return true)
    end

    return false
end

"""
    resetEngine!(p4MassConstantContinuationEngine, solution1, solution2)

Return P4 mass continuation engine with reset data

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
- `solution1::BCR4BP12MultipleShooterProblem`: First member of family
- `solution2::BCR4BP12MultipleShooterProblem`: Second member of family
"""
function resetEngine!(p4MassContinuationEngine::P4MassContinuationEngine, solution1::MBD.BCR4BP12MultipleShooterProblem, solution2::MBD.BCR4BP12MultipleShooterProblem)
    p4MassContinuationEngine.dataInProgress = MBD.BCR4BP12ContinuationData(solution1, solution2)
end

"""
    tryConverging!(p4MassContinuationEngine)

Return updated P4 mass continuation engine object

# Arguments
- `p4MassContinuationEngine::P4MassContinuationEngine`: P4 mass continuation engine object
"""
function tryConverging!(p4MassContinuationEngine::P4MassContinuationEngine)
    updateStepSize!(p4MassContinuationEngine.stepSizeGenerator, p4MassContinuationEngine.dataInProgress)
    p4MassContinuationEngine.printProgress && println("Current step size: $(p4MassContinuationEngine.dataInProgress.currentStepSize)")
    p4MassContinuationEngine.dataInProgress.nextGuess = deepClone(p4MassContinuationEngine.dataInProgress.previousSolution)
    setFreeVariableVector!(p4MassContinuationEngine.dataInProgress.nextGuess, getFreeVariableVector!(p4MassContinuationEngine.dataInProgress.previousSolution)+p4MassContinuationEngine.dataInProgress.fullStep.*p4MassContinuationEngine.dataInProgress.currentStepSize)
    constrainNextGuess!(p4MassContinuationEngine, p4MassContinuationEngine.dataInProgress)
    try
        twoPreviousConvergedSolution::MBD.BCR4BP12MultipleShooterProblem = p4MassContinuationEngine.dataInProgress.twoPreviousSolution
        previousConvergedSolution::MBD.BCR4BP12MultipleShooterProblem = p4MassContinuationEngine.dataInProgress.previousSolution
        p4MassContinuationEngine.dataInProgress.twoPreviousSolution = deepClone(p4MassContinuationEngine.dataInProgress.previousSolution)
        p4MassContinuationEngine.dataInProgress.previousSolution = solve!(p4MassContinuationEngine.corrector, p4MassContinuationEngine.dataInProgress.nextGuess)
        p4MassContinuationEngine.dataInProgress.converging = true
        if abs(LinearAlgebra.norm(getFreeVariableVector!(p4MassContinuationEngine.dataInProgress.previousSolution))-LinearAlgebra.norm(getFreeVariableVector!(p4MassContinuationEngine.dataInProgress.twoPreviousSolution))) > abs(p4MassContinuationEngine.dataInProgress.currentStepSize)*50
            println("Entered if")
            p4MassContinuationEngine.dataInProgress.converging = false
            p4MassContinuationEngine.printProgress && println("Solution outside of trust region: delta = $(abs(LinearAlgebra.norm(getFreeVariableVector!(p4MassContinuationEngine.dataInProgress.previousSolution))-LinearAlgebra.norm(getFreeVariableVector!(p4MassContinuationEngine.dataInProgress.twoPreviousSolution))))")
        else
            println("Entered else")
            for jumpCheck::MBD.AbstractContinuationJumpCheck in p4MassContinuationEngine.jumpChecks
                if typeof(jumpCheck) == MBD.BoundingBoxJumpCheck
                    println("Check reached")
                    for (index::MBD.Variable, value::Int16) in p4MassContinuationEngine.dataInProgress.previousSolution.freeVariableIndexMap
                        if index.name == jumpCheck.paramName
                            addBounds!(jumpCheck, p4MassContinuationEngine.dataInProgress.previousSolution, index, jumpCheck.paramBounds)
                            p4MassContinuationEngine.dataInProgress.converging = isFamilyMember(jumpCheck, p4MassContinuationEngine.dataInProgress)
                            (!p4MassContinuationEngine.dataInProgress.converging && p4MassContinuationEngine.printProgress) && println("Solution jumped")
                            removeBounds!(jumpCheck, p4MassContinuationEngine.dataInProgress.previousSolution, index)
                            break
                        end
                    end
                end
                !p4MassContinuationEngine.dataInProgress.converging && break
            end
        end
        if p4MassContinuationEngine.dataInProgress.converging
            p4MassContinuationEngine.dataInProgress.numIterations = p4MassContinuationEngine.corrector.recentIterationCount
        else
            p4MassContinuationEngine.dataInProgress.twoPreviousSolution = twoPreviousConvergedSolution
            p4MassContinuationEngine.dataInProgress.previousSolution = previousConvergedSolution
        end
    catch
        p4MassContinuationEngine.dataInProgress.converging = false
        p4MassContinuationEngine.printProgress && println("Failed to converge")
    end
end
