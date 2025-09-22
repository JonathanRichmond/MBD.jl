"""
Keplerian system data wrapper

Author: Jonathan Richmond
C: 9/17/25
"""

import MBD: KSystemData

export getMassParameter, getNumPrimaries

# """
#     getCharLength(systemData)

# Return CR3BP characteristic length

# # Arguments
# - `systemData::CR3BPSystemData`: CR3BP system data object
# """
# function getCharLength(systemData::CR3BPSystemData)
#     return systemData.primaryData[2].orbitRadius
# end

# """
#     getCharMass(systemData)

# Return CR3BP characteristic mass

# # Arguments
# - `systemData::CR3BPSystemData`: CR3BP system data object
# """
# function getCharMass(systemData::CR3BPSystemData)
#     totalGravParam::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam
    
#     return totalGravParam/GRAVITY
# end

# """
#     getCharTime(systemData)

# Return CR3BP characteristic time

# # Arguments
# - `systemData::CR3BPSystemData`: CR3BP system data object
# """
# function getCharTime(systemData::CR3BPSystemData)
#     lstar::Float64 = getCharLength(systemData)
#     totalGravParam::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam
    
#     return sqrt(lstar^3/totalGravParam)
# end

"""
    getMassParameter(systemData)

Return Keplerian system mass parameter

# Arguments
- `systemData::KSystemData`: Keplerian system data object
"""
function getMassParameter(systemData::KSystemData)
    return systemData.primaryData.gravParam
end

"""
    getNumPrimaries(systemData)

Return number of primaries

# Arguments
- `systemData::KSystemData`: Keplerian system data object
"""
function getNumPrimaries(systemData::KSystemData)
    return Int16(length(systemData.primaryData))
end

"""
    shallowClone(systemData)

Return copy of Keplerian system data object

# Arguments
- `systemData::KSystemData`: Keplerian system data object
"""
function shallowClone(systemData::KSystemData)
    return KSystemData(systemData.primaryName)
end
