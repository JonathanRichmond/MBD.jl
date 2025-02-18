"""
BCR4BP system data wrapper

Author: Jonathan Richmond
C: 2/18/25
"""

import MBD: BCR4BPSystemData, GRAVITY

export getNumPrimaries, get12CharLength, get12CharMass, get12CharTime, get12MassRatio, get4Distance
export get4Mass, get41CharLength, get41CharMass, get41CharTime, get41MassRatio

"""
    getNumPrimaries(systemData)

Return number of primaries

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function getNumPrimaries(systemData::BCR4BPSystemData)
    return Int16(length(systemData.primaryData)-1)
end

"""
    get12CharLength(systemData)

Return P1-P2 CR3BP characteristic length

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get12CharLength(systemData::BCR4BPSystemData)
    return systemData.primaryData[2].orbitRadius
end

"""
    get12CharMass(systemData)

Return P1-P2 CR3BP characteristic mass

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get12CharMass(systemData::BCR4BPSystemData)
    gravParam12::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam
    
    return gravParam12/GRAVITY
end

"""
    get12CharTime(systemData)

Return P1-P2 CR3BP characteristic time

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get12CharTime(systemData::BCR4BPSystemData)
    lstar12::Float64 = get12CharLength(systemData)
    gravParam12::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam
    
    return sqrt(lstar12^3/gravParam12)
end

"""
    get12MassRatio(systemData)

Return P1-P2 CR3BP system mass ratio

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get12MassRatio(systemData::BCR4BPSystemData)
    gravParam12::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam

    return systemData.primaryData[2].gravParam/gravParam12
end

"""
    get4Distance(systemData)

Return P4 distance from B1 [ndim]

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get4Distance(systemData::BCR4BPSystemData)
    return systemData.primaryData[4].orbitRadius/get12CharLength(systemData)
end

"""
    get4Mass(systemData)

Return P4 mass [ndim]

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get4Mass(systemData::BCR4BPSystemData)
    return systemData.primaryData[3].gravParam/GRAVITY/get12CharMass(systemData)
end

"""
    get41CharLength(systemData)

Return P4-B1 CR3BP characteristic length

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get41CharLength(systemData::BCR4BPSystemData)
    return get12CharLength(systemData)*get4Distance(systemData)
end

"""
    get41CharMass(systemData)

Return P4-B1 CR3BP characteristic mass

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get41CharMass(systemData::BCR4BPSystemData) 
    return get12CharMass(systemData)*(get4Mass(systemData)+1)
end

"""
    get41CharTime(systemData)

Return P4-B1 CR3BP characteristic time

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get41CharTime(systemData::BCR4BPSystemData)
    lstar41::Float64 = get41CharLength(systemData)
    gravParam41::Float64 = systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam+systemData.primaryData[3].gravParam
    
    return sqrt(lstar41^3/gravParam41)
end

"""
    get41MassRatio(systemData)

Return P4-B1 CR3BP system mass ratio

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
function get41MassRatio(systemData::BCR4BPSystemData)
    return get12CharMass(systemData)/get41CharMass(systemData)
end
