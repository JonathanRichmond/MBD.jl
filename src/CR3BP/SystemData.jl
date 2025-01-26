"""
CR3BP system data wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 1/13/25
"""

import MBD: CR3BPSystemData, GRAVITY

export getCharLength, getCharMass, getCharTime, getMassRatio, getNumPrimaries

"""
    getCharLength(systemData)

Return CR3BP characteristic length

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system data object
"""
function getCharLength(systemData::CR3BPSystemData)
    return systemData.primaryData[2].orbitRadius
end

"""
    getCharMass(systemData)

Return CR3BP characteristic mass

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system data object
"""
function getCharMass(systemData::CR3BPSystemData)
    totalGravParam::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam
    
    return totalGravParam/GRAVITY
end

"""
    getCharTime(systemData)

Return CR3BP characteristic time

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system data object
"""
function getCharTime(systemData::CR3BPSystemData)
    lstar::Float64 = getCharLength(systemData)
    totalGravParam::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam
    
    return sqrt(lstar^3/totalGravParam)
end

"""
    getMassRatio(systemData)

Return CR3BP system mass ratio

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system data object
"""
function getMassRatio(systemData::CR3BPSystemData)
    totalGravParam::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam

    return systemData.primaryData[2].gravParam/totalGravParam
end

"""
    getNumPrimaries(systemData)

Return number of primaries

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system data object
"""
function getNumPrimaries(systemData::CR3BPSystemData)
    return Int16(length(systemData.primaryData))
end
