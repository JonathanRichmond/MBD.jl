"""
Number of steps continuation end check wrapper

Author: Jonathan Richmond
C: 1/5/23
U: 6/23/25
"""

import Logging
import MBD: NumberStepsContinuationEndCheck

export isContinuationDone

"""
    isContinuationDone(numStepsCheck, data)

Return true if continuation is done

# Arguments
- `numStepsCheck::NumberStepsContinuationEndCheck`: Number of steps continuation end check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isContinuationDone(numStepsCheck::NumberStepsContinuationEndCheck, data::MBD.CR3BPContinuationData)
    numSteps::Int64 = getNumSteps(data)
    maxSteps::Int64 = numStepsCheck.maxSteps

    Logging.@debug "Checking if number of continuation steps is reached"

    if numSteps >= maxSteps
        Logging.@info "Number of continuation steps reached: $numSteps ≥ $maxSteps"
        println("Number of continuation steps reached!")

        return true
    end

    Logging.@debug "Continuation ongoing: $numSteps / $maxSteps steps completed"
    return false
end

"""
    isContinuationDone(numStepsCheck, data)

Return true if continuation is done

# Arguments
- `numStepsCheck::NumberStepsContinuationEndCheck`: Number of steps continuation end check object
- `data::BCR4BP12ContinuationData`: BCR4BP P1-P2 continuation data object
"""
function isContinuationDone(numStepsCheck::NumberStepsContinuationEndCheck, data::MBD.BCR4BP12ContinuationData)
    numSteps::Int64 = getNumSteps(data)
    maxSteps::Int64 = numStepsCheck.maxSteps

    Logging.@debug "Checking if number of continuation steps is reached"

    if numSteps >= maxSteps
        Logging.@info "Number of continuation steps reached: $numSteps ≥ $maxSteps"
        println("Number of continuation steps reached!")

        return true
    end

    Logging.@debug "Continuation ongoing: $numSteps / $maxSteps steps completed"
    return false
end
