"""
SPICE functions

Author: Jonathan Richmond
C: 9/14/23
"""

import SPICE

export getEphemerides

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
    SPICE.furnsh("MBD.jl/src/spice/kernels/naif0012.tls", "MBD.jl/src/spice/kernels/de440.bsp", "MBD.jl/src/spice/kernels/mar097.bsp")
    epoch::Float64 = SPICE.str2et(initialEpoch)
    epochs::Vector{Float64} = epoch.+times
    ephemerisStates::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(epochs))
    ephemerisTimes::Vector{Float64} = Vector{Float64}(undef, length(epochs))
    for e::Int64 in 1:length(epochs)
        (state::Vector{Float64}, time::Float64) = SPICE.spkezr(targetBody, epochs[e], "ECLIPJ2000", "NONE", "SUN")
        ephemerisStates[e] = state
        ephemerisTimes[e] = time
    end
    SPICE.kclear()
    
    return (ephemerisStates, ephemerisTimes)
end
