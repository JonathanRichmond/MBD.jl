"""
CR3BP continuation family wrapper

Author: Jonathan Richmond
C: 1/10/25
"""

import MBD: CR3BPContinuationFamily

export getNumMembers

"""
    getNumMembers(continuationFamily)

Return number of continuation family members

# Arguments
- `continuationFamily::CR3BPContinuationFamily`: CR3BP continuation family object
"""
function getNumMembers(continuationFamily::CR3BPContinuationFamily)
    return length(continuationFamily.segments)
end
