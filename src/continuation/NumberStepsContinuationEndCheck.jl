"""
Number of steps continuation end check wrapper

Author: Jonathan Richmond
C: 1/5/23
U: 9/8/23
"""

import MBD: NumberStepsContinuationEndCheck

export isContinuationDone

"""
    isContinuationDone(numberStepsContinuationEndCheck, data)

Return true if continuation is done

# Arguments
- `numberStepsContinuationEndCheck::NumberStepsContinuationEndCheck`: Number of steps continuation end check object
- `data::ContinuationData`: Continuation data object
"""
function isContinuationDone(numberStepsContinuationEndCheck::NumberStepsContinuationEndCheck, data::MBD.ContinuationData)
    (data.stepCount >= numberStepsContinuationEndCheck.maxSteps) && println("Number of continuation steps reached!")

    return (data.stepCount >= numberStepsContinuationEndCheck.maxSteps)
end
