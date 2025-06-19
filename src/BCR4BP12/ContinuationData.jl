"""
BCR4BP P1-P2 continuation data wrapper

Author: Jonathan Richmond
C: 6/18/25
"""

import MBD: BCR4BP12ContinuationData

export getNumSteps

"""
    getNumSteps(continuationData)

Return number of continuation steps

# Arguments
- `continuationData::BCR4BP12ContinuationData`: BCR4BP P1-P2 continuation data object
"""
function getNumSteps(continuationData::BCR4BP12ContinuationData)
    return getNumMembers(continuationData.family)
end
