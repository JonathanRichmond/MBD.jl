"""
BodyName wrapper

Author: Jonathan Richmond
C: 9/1/22
U: 7/29/23
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
    return SPICE.bods2c(bodyName.name)
end