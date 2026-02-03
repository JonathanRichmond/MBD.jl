"""
SPICE functions

Author: Jonathan Richmond
C: 9/14/23
U: 2/3/26
"""

import Ephemerides, SPICE

export getEphemerides, getSunEphemerides

"""
    getEphemerides(initialEpoch, times, targetBody, observerBody, frame)

Return ephemerides

# Arguments
- `initialEpoch::String`: Initial epoch
- `times::Vector{Float64}`: Times since initial epoch
- `targetBody::String`: Target body for ephemerides
- `observerBody::String`: Body for ephemerides reference
- `frame::String`: Reference frame
"""
function getEphemerides(initialEpoch::String, times::Vector{Float64}, targetBody::String, observerBody::String, frame::String)
    epoch::Float64 = SPICE.str2et(initialEpoch)
    epochs::Vector{Float64} = epoch.+times
    numStates::Int16 = length(epochs)
    ephemerisStates::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numStates)
    ephemerisTimes::Vector{Float64} = Vector{Float64}(undef, numStates)
    for e::Int64 in 1:length(epochs)
        (state::Vector{Float64}, time::Float64) = SPICE.spkezr(targetBody, epochs[e], frame, "NONE", observerBody)
        ephemerisStates[e] = state
        ephemerisTimes[e] = time
    end
    
    return (ephemerisStates, ephemerisTimes)
end

"""
    getEphemerides(eph, initialEpochTime, times, targetBodyID, observerBodyID, referenceBodyID; frame)

Return ephemerides

# Arguments
- `eph::EphemerisProvider`: Ephemeris provider object
- `initialEpochTime::Float64`: Initial epoch time [s]
- `times::Vector{Float64}`: Times since initial epoch
- `targetBodyID::Int16`: Target body SPICE ID for ephemerides
- `observerBodyID::Int16`: Observer body SPICE ID for ephemerides
- `referenceBodyID::Int16`: Body SPICE ID for ephemerides reference (common to both)
- `frame::String`: Reference frame (default = "ECLIPJ2000")
"""
function getEphemerides(eph::Ephemerides.EphemerisProvider, initialEpochTime::Float64, times::Vector{Float64}, targetBodyID::Int16, observerBodyID::Int16, referenceBodyID::Int64; frame::String = "ECLIPJ2000")
    ephemerisTimes::Vector{Float64} = initialEpochTime .+ times
    if frame == "ECLIPJ2000"
        i::Float64 = 23.43929111*pi/180
        R::Matrix{Float64} = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]
        N::Matrix{Float64} = [R zeros(Float64, (3,3)); zeros(Float64, (3,3)) R]
        ephemerisStates::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(ephemerisTimes))
        for et::Int64 in 1:length(ephemerisTimes)
            eqState::Vector{Float64} = Ephemerides.ephem_vector6(eph, referenceBodyID, Int64(targetBodyID), ephemerisTimes[et])-Ephemerides.ephem_vector6(eph, referenceBodyID, Int64(observerBodyID), ephemerisTimes[et])
            ephemerisStates[et] = N*eqState
        end
    elseif frame == "J2000"
        ephemerisStates = [Ephemerides.ephem_vector6(eph, referenceBodyID, targetBodyID, et)-Ephemerides.ephem_vector6(eph, referenceBodyID, observerBodyID, et) for et in ephemerisTimes]
    else
        throw(ArgumentError("Frame not supported"))
    end

    return (ephemerisStates, ephemerisTimes)
end

"""
    getSunEphemerides(eph, initialEpochTime, times, targetBodyID, observerBodyID; frame)

Return ephemerides relative to the Sun

# Arguments
- `eph::EphemerisProvider`: Ephemeris provider object
- `initialEpochTime::Float64`: Initial epoch time [s]
- `times::Vector{Float64}`: Times since initial epoch
- `targetBodyID::Int16`: Target body SPICE ID for ephemerides
- `observerBodyID::Int16`: Observer body SPICE ID for ephemerides
- `frame::String`: Reference frame (default = "ECLIPJ2000")
"""
function getSunEphemerides(eph::Ephemerides.EphemerisProvider, initialEpochTime::Float64, times::Vector{Float64}, targetBodyID::Int16, observerBodyID::Int16; frame::String = "ECLIPJ2000")
    ephemerisTimes::Vector{Float64} = initialEpochTime .+ times
    ephemerisStates::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(ephemerisTimes))
    for et::Int64 in 1:length(ephemerisTimes)
        bodyState::Vector{Float64} = Ephemerides.ephem_vector6(eph, firstdigit(targetBodyID), Int64(targetBodyID), ephemerisTimes[et])
        barycenterState::Vector{Float64} = Ephemerides.ephem_vector6(eph, 0, firstdigit(targetBodyID), ephemerisTimes[et])
        SunState::Vector{Float64} = Ephemerides.ephem_vector6(eph, observerBodyID, 0, ephemerisTimes[et])
        eqState::Vector{Float64} = bodyState-barycenterState-SunState
        if frame == "ECLIPJ2000"
            i::Float64 = 23.43929111*pi/180
            R::Matrix{Float64} = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]
            N::Matrix{Float64} = [R zeros(Float64, (3,3)); zeros(Float64, (3,3)) R]
            ephemerisStates[et] = N*eqState
        elseif frame == "J2000"
            ephemerisStates[et] = eqState
        else
            throw(ArgumentError("Frame not supported"))
        end
    end
    
    return (ephemerisStates, ephemerisTimes)
end

