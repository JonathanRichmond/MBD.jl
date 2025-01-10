"""
CR3BP continuation data wrapper

Author: Jonathan Richmond
C: 1/10/25
"""

import MBD: CR3BPContinuationData

export getNumSteps

"""
    getNumSteps(continuationData)

Return number of continuation steps

# Arguments
- `continuationData::CR3BPContinuationData`: CR3BP continuation data object
"""
function getNumSteps(continuationData::CR3BPContinuationData)
    return getNumMembers(continuationData.family)
end
