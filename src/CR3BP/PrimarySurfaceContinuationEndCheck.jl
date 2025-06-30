"""
CR3BP primary surface continuation end check wrapper

Author: Jonathan Richmond
C: 6/30/25
"""

import LinearAlgebra, Logging
import MBD: CR3BPPrimarySurfaceContinuationEndCheck

export isContinuationDone

"""
    isContinuationDone(boundsCheck, data)

Return true if continuation is done

# Arguments
- `boundsCheck::CR3BPPrimarySurfaceContinuationEndCheck`: CR3BP primary surface continuation end check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isContinuationDone(primarySurfaceCheck::CR3BPPrimarySurfaceContinuationEndCheck, data::MBD.CR3BPContinuationData)
    primaryPos::Vector{Float64} = getPrimaryState(primarySurfaceCheck.dynamicsModel, primarySurfaceCheck.primary)[1:3]
    primaryRad::Float64 = primarySurfaceCheck.dynamicsModel.systemData.primaryData[primarySurfaceCheck.primary].bodyRadius/getCharLength(primarySurfaceCheck.dynamicsModel)
    Logging.@debug "Checking if primary surface is reached"

    d::Float64 = LinearAlgebra.norm(data.previousSolution.nodes[1].state.data[1:3]-primaryPos)
    if d <= primaryRad
        Logging.@info "Primary surface reached: $d ≤ $primaryRad"
        println("Primary surface reached!")

        return true
    end

    Logging.@debug "Continuation ongoing: primary surface not yet reached"
    return false
end
