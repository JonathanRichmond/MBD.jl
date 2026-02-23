"""
CR3BP manifold wrapper

Author: Jonathan Richmond
C: 1/25/25
U: 2/23/26
"""

import DifferentialEquations
import MBD: CR3BPManifold

export getJacobiConstant, stopCrashes

"""
    getJacobiConstant(manifold)

Return Jacobi constant

# Arguments
- `manifold::CR3BPManifold`: CR3BP manifold object
"""
function getJacobiConstant(manifold::CR3BPManifold)
    return getJacobiConstant(manifold.periodicOrbit)
end

"""
    stopCrashes(manifold)

Return manifold arcs, stopping propagation when a primary is encountered

# Arguments
- `manifold::CR3BPManifold`: CR3BP manifold object
"""
function stopCrashes(manifold::CR3BPManifold)
    propagator = MBD.Propagator()
    crashEvent = DifferentialEquations.VectorContinuousCallback(primaryDistanceCondition2, terminateAffectIndex!, 2)
    P1Radius::Float64 = manifold.periodicOrbit.dynamicsModel.systemData.primaryData[1].bodyRadius/getCharLength(manifold.periodicOrbit.dynamicsModel)
    P2Radius::Float64 = manifold.periodicOrbit.dynamicsModel.systemData.primaryData[2].bodyRadius/getCharLength(manifold.periodicOrbit.dynamicsModel)
    manifoldArcs::Vector{MBD.CR3BPManifoldArc} = Vector{MBD.CR3BPManifoldArc}(undef, length(manifold.initialConditions))
    for a::Int64 = 1:length(manifold.initialConditions)
        arc::MBD.CR3BPArc = propagateWithEvent(propagator, crashEvent, real(manifold.initialConditions[a]), [0.0, manifold.TOF], manifold.periodicOrbit.dynamicsModel, [manifold.periodicOrbit.dynamicsModel, P1Radius, P2Radius])
        manifoldArcs[a] = MBD.CR3BPManifoldArc(manifold.periodicOrbit, manifold.orbitTimes[a], manifold.direction, manifold.ds[a], manifold.initialConditions[a], getTimeByIndex(arc, -1))
    end

    return manifoldArcs
end
