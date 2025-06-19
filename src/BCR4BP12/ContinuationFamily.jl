"""
BCR5BP P1-P2 continuation family wrapper

Author: Jonathan Richmond
C: 6/18/25
"""

import MBD: BCR4BP12ContinuationFamily

export getNumMembers

"""
    getNumMembers(continuationFamily)

Return number of continuation family members

# Arguments
- `continuationFamily::BCR4BP12ContinuationFamily`: BCR4BP P1-P2 continuation family object
"""
function getNumMembers(continuationFamily::BCR4BP12ContinuationFamily)
    return length(continuationFamily.segments)
end
