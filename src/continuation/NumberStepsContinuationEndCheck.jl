"""
Number of steps continuation end check wrapper

Author: Jonathan Richmond
C: 1/5/23
U: 1/26/25
"""

import MBD: NumberStepsContinuationEndCheck

export isContinuationDone

"""
    isContinuationDone(numberStepsContinuationEndCheck, data)

Return true if continuation is done

# Arguments
- `numberStepsContinuationEndCheck::NumberStepsContinuationEndCheck`: Number of steps continuation end check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isContinuationDone(numberStepsContinuationEndCheck::NumberStepsContinuationEndCheck, data::MBD.CR3BPContinuationData)
    (getNumSteps(data) >= numberStepsContinuationEndCheck.maxSteps) && println("Number of continuation steps reached!")

    return (getNumSteps(data) >= numberStepsContinuationEndCheck.maxSteps)
end
