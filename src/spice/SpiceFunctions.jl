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
    SPICE.furnsh("MBD.jl/src/spice/kernels/naif0012.tls")#, "kernels/de440.bsp", "kernels/mar097.bsp")
    epoch = SPICE.str2et(initialEpoch)
    epochs = epoch+times
    (states, times) = SPICE.spkezr(targetBody, epochs, frame, "NONE", observerBody)
    SPICE.kclear()
    
    return (states, times)
end
