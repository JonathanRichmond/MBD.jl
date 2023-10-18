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
    epoch::Float64 = SPICE.str2et(initialEpoch)
    epochs::Vector{Float64} = epoch.+times
    ephemerisStates::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(epochs))
    for e::Int64 in 1:length(epochs)
        (state::Vector{Float64}, time::Float64) = SPICE.spkezr(targetBody, epochs[e], frame, "NONE", observerBody)
        ephemerisStates[e] = state
    end
    
    return ephemerisStates
end
