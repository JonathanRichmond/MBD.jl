"""
Body name wrapper

Author: Jonathan Richmond
C: 9/1/22
U: 1/9/25
"""

import SPICE
import MBD: BodyName

export getIDCode

"""
    getIDCode(bodyName)

Return body ID code

# Arguments
- `bodyName::BodyName`: BodyName object
"""
function getIDCode(bodyName::BodyName)
    return Int16(SPICE.bods2c(bodyName.name))
end
