"""
CR3BP multiple shooter continuation engine wrapper

Author: Jonathan Richmond
C: 4/1/25
"""

import StaticArrays
import MBD: CR3BPMultipleShooterContinuationEngine

export addEndCheck!, addJumpCheck!, computeStateStep, computeTimeStep, convergeInitialSolution
export endContinuation, resetEngine!

"""
    addEndCheck!(multipleShooterContinuationEngine, endCheck)

Return CR3BP multiple shooter continuation engine with updated end checks

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `endCheck::AbstractContinuationEndCheck`: Continuation end check
"""
function addEndCheck!(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, endCheck::MBD.AbstractContinuationEndCheck)
    push!(multipleShooterContinuationEngine.endChecks, endCheck)
end

"""
    addJumpCheck!(multipleShooterContinuationEngine, jumpCheck)

Return CR3BP multiple shooter continuation engine with updated jump checks

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `jumpCheck::AbstractContinuationJumpCheck`: Continuation jump check
"""
function addJumpCheck!(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, jumpCheck::MBD.AbstractContinuationJumpCheck)
    push!(multipleShooterContinuationEngine.jumpChecks, jumpCheck)
end

"""
    computeStateStep(multipleShooterContinuationEngine, data)

Return update step for initial state

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function computeStateStep(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, data::MBD.CR3BPContinuationData)
    state1::Vector{Float64} = data.twoPreviousSolution.nodes[1].state.data[1:6]
    state2::Vector{Float64} = data.previousSolution.nodes[1].state.data[1:6]

    return (state2-state1)./data.currentStepSize
end

"""
    computeTimeStep(multipleShooterContinuationEngine, data)

Return update step for final time

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function computeTimeStep(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, data::MBD.CR3BPContinuationData)
    time1::Float64 = 0.0
    time2::Float64 = 0.0
    for s = 1:length(data.twoPreviousSolution.segments)
        time1 += data.twoPreviousSolution.segments[s].TOF.data[1]
        time2 += data.previousSolution.segments[s].TOF.data[1]
    end
    time1 *= 2
    time2 *= 2

    return (time2-time1)/data.currentStepSize
end

"""
    convergeInitialSolution(multipleShooterContinuationEngine, initialGuess)

Return converged initial solution

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `initialGuess::CR3BPMultipleShooterProblem`: Initial family member
"""
function convergeInitialSolution(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, initialGuess::MBD.CR3BPMultipleShooterProblem)
    return solve!(multipleShooterContinuationEngine.corrector, initialGuess)
end

"""
    endContinuation(multipleShooterContinuationEngine, data)

Return true if continuation should end

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`:CR3BP multiple shooter continuation engine object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function endContinuation(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, data::MBD.CR3BPContinuationData)
    for endCheck::MBD.AbstractContinuationEndCheck in multipleShooterContinuationEngine.endChecks
        isContinuationDone(endCheck, data) && (return true)
    end

    return false
end

"""
    resetEngine!(multipleShooterContinuationEngine, solution1, solution2)

Return CR3BP multiple shooter continuation engine with reset data

# Arguments
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `solution1::CR3BPMultipleShooterProblem`: First member of family
- `solution2::CR3BPMultipleShooterProblem`: Second member of family
"""
function resetEngine!(multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine, solution1::MBD.CR3BPMultipleShooterProblem, solution2::MBD.CR3BPMultipleShooterProblem)
    multipleShooterContinuationEngine.dataInProgress = MBD.CR3BPContinuationData(solution1, solution2)
end
