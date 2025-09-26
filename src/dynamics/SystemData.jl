"""
System data wrapper

Author: Jonathan Richmond
C: 7/24/25
"""

import Logging
import MBD: SystemData

export getCharLengths, getCharMasses, getCharTimes, getMassParams, getNumPrimaries, getSupParams


"""
    getCharLengths(systemData::SystemData) -> Vector{Float64}

Return system characteristic length(s)

# Arguments
- `systemData::SystemData`: SystemData object

# Outputs
- `Vector{Float64}`: System characteristic length(s)
    - For CR3BP: `[lstar]`
    - For BCR4BP: `[lstar12, lstar41]`
"""
function getCharLengths(systemData::SystemData)
    numPrimaries::Int64 = getNumPrimaries(systemData)
    if numPrimaries < 2
        err1::String = "System must have at least two primaries to define characteristic length(s)"
        Logging.@error err1
        throw(ArgumentError(err1))
    end

    Logging.@debug "Retrieving $(systemData.modelType) system characteristic length(s)"
 
    if systemData.modelType == MBD.CR3BP
        lstar::Float64 = systemData.primaryData[2].orbitRadius
        Logging.@debug "Returning one characteristic length"

        return [lstar]
    elseif systemData.modelType == MBD.BCR4BP
        lstar12::Float64 = systemData.primaryData[3].orbitRadius
        supParams::Vector{Float64} = getSupParams(systemData)
        if length(supParams) != 2
            err2::String = "Insufficient supplemental parameters for BCR4BP"
            Logging.@error err2
            throw(ErrorException(err2))
        end
        lstar41::Float64 = lstar12*supParams[1]
        Logging.@debug "Returning two characteristic lengths"

        return [lstar12, lstar41]
    else
        err3::String = "Unsupported number of primaries: $numPrimaries; expected 2 or 3"
        Logging.@error err3
        throw(ErrorException(err3))
    end
end

"""
    getCharMasses(systemData::SystemData) -> Vector{Float64}

Return system characteristic mass(es)

# Arguments
- `systemData::SystemData`: SystemData object

# Outputs
- `Vector{Float64}`: System characteristic mass(es)
    - For CR3BP: `[mstar]`
    - For BCR4BP: `[mstar12, mstar41]`
"""
function getCharMasses(systemData::SystemData)
    numPrimaries::Int64 = getNumPrimaries(systemData)
    if numPrimaries < 2
        err1::String = "System must have at least two primaries to define characteristic mass(es)"
        Logging.@error err1
        throw(ArgumentError(err1))
    end

    Logging.@debug "Retrieving $(systemData.modelType) system characteristic mass(es)"
 
    if systemData.modelType == MBD.CR3BP
        mstar::Float64 = (systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam)/GRAVITY
        Logging.@debug "Returning one characteristic mass"

        return [mstar]
    elseif systemData.modelType == MBD.BCR4BP
        mstar12::Float64 = (systemData.primaryData[2].gravParam+systemData.primaryData[3].gravParam)/GRAVITY
        mstar41::Float64 = mstar12+systemData.primaryData[1].gravParam/GRAVITY
        Logging.@debug "Returning two characteristic masses"

        return [mstar12, mstar41]
    else
        err2::String = "Unsupported number of primaries: $numPrimaries; expected 2 or 3"
        Logging.@error err2
        throw(ErrorException(err2))
    end
end

"""
    getCharTimes(systemData::SystemData) -> Vector{Float64}

Return system characteristic time(s)

# Arguments
- `systemData::SystemData`: SystemData object

# Outputs
- `Vector{Float64}`: System characteristic time(s)
    - For CR3BP: `[tstar]`
    - For BCR4BP: `[tstar12, tstar41]`
"""
function getCharTimes(systemData::SystemData)
    numPrimaries::Int64 = getNumPrimaries(systemData)
    if numPrimaries < 2
        err1::String = "System must have at least two primaries to define characteristic time(s)"
        Logging.@error err1
        throw(ArgumentError(err1))
    end

    Logging.@debug "Retrieving $(systemData.modelType) system characteristic time(s)"
 
    if systemData.modelType == MBD.CR3BP
        tstar::Float64 = sqrt(getCharLengths(systemData)[1]^3/(systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam))
        Logging.@debug "Returning one characteristic time"

        return [tstar]
    elseif systemData.modelType == MBD.BCR4BP
        tstar12::Float64 = sqrt(getCharLengths(systemData)[1]^3/(systemData.primaryData[2].gravParam+systemData.primaryData[3].gravParam))
        tstar41::Float64 = sqrt(getCharLengths(systemData)[2]^3/(systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam+systemData.primaryData[3].gravParam))
        Logging.@debug "Returning two characteristic times"

        return [tstar12, tstar41]
    else
        err2::String = "Unsupported number of primaries: $numPrimaries; expected 2 or 3"
        Logging.@error err2
        throw(ErrorException(err2))
    end
end

"""
    getMassParams(systemData::SystemData) -> Vector{Float64}

Return system mass parameter(s)

# Arguments
- `systemData::SystemData`: SystemData object

# Outputs
- `Vector{Float64}`: System mass parameter(s)
    - For CR3BP: `[mu]`
    - For BCR4BP: `[mu12, mu41]`
"""
function getMassParams(systemData::SystemData)
    numPrimaries::Int64 = getNumPrimaries(systemData)
    if numPrimaries < 1
        err1::String = "System must have at least one primary to define mass parameter(s)"
        Logging.@error err1
        throw(ArgumentError(err1))
    end

    Logging.@debug "Retrieving $(systemData.modelType) system mass parameter(s)"
 
    if systemData.modelType == MBD.CR3BP
        mu::Float64 = systemData.primaryData[2].gravParam/(systemData.primaryData[1].gravParam+systemData.primaryData[2].gravParam)
        Logging.@debug "Returning one mass ratio"

        return [mu]
    elseif systemData.modelType == MBD.BCR4BP
        mu12::Float64 = systemData.primaryData[3].gravParam/(systemData.primaryData[2].gravParam+systemData.primaryData[3].gravParam)
        mstars::Vector{Float64} = getCharMasses(systemData)
        mu41::Float64 = mstars[1]/mstars[2]
        Logging.@debug "Returning two mass ratios"

        return [mu12, mu41]
    else
        err2::String = "Unsupported number of primaries: $numPrimaries; expected 2 or 3"
        Logging.@error err2
        throw(ErrorException(err2))
    end
end

"""
    getNumPrimaries(systemData::SystemData) -> Int64

Return number of system primaries

# Arguments
- `systemData::SystemData`: SystemData object

# Outputs
- `Int64`: Number of system primaries
"""
function getNumPrimaries(systemData::SystemData)
    Logging.@debug "Returning number of primaries in $(systemData.modelType) system"

    return length(systemData.primaryData)
end

"""
    getSupParams(systemData::SystemData) -> Vector{Float64}

Return any supplemental system parameters

# Arguments
- `systemData::SystemData`: SystemData object

# Outputs
- `Vector{Float64}`: Supplemental system parameters
    - For BCR4BP: `[a4, m4]`
"""
function getSupParams(systemData::SystemData)
    Logging.@debug "Retrieving $(systemData.modelType) system supplemental parameters"

    if systemData.modelType == MBD.BCR4BP
        a4::Float64 = systemData.primaryData[2].orbitRadius/systemData.primaryData[3].orbitRadius
        m4::Float64 = systemData.primaryData[1].gravParam/GRAVITY/getCharMasses(systemData)[1]
        Logging.@debug "Returning nondimensionalized P4 distance and mass"

        return [a4, m4]
    else
        err1::String = "$(systemData.modelType) system has no supplemental parameters"
        Logging.@error err1
        throw(ErrorException(err1))
    end
end
