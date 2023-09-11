"""
CR3BP system data wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 7/30/23
"""

import MBD: CR3BPSystemData

export getMassRatio

"""
    getMassRatio(systemData)

Return CR3BP system mass ratio

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system object
"""
function getMassRatio(systemData::CR3BPSystemData)
    return systemData.params[1]
end
