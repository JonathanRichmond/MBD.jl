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
    ephemerisStates::Vector{Float64} = Vector{Float64}(undef, length(times))
    ephemerisTimes::Vector{Float64} = Vector{Float64}(undef, length(times))
    for e::Int64 in length(times)
        (state::Vector{Float64}, time::Float64) = SPICE.spkezr(targetBody, epochs[e], frame, "NONE", observerBody)
        ephemerisStates[e] = state
        ephemerisTimes[e] = time
    end
    SPICE.kclear()
    
    return (states, times)
end
