"""
BCR4BP P1-P2 manifold wrapper

Author: Jonathan Richmond
C: 6/16/25
"""

import DifferentialEquations
import MBD: BCR4BP12Manifold

export stopCrashes

"""
    stopCrashes(manifold)

Return manifold arcs, stopping propagation when a primary is encountered

# Arguments
- `manifold::BCR4BP12Manifold`: BCR4BP P1-P2 manifold object
"""
function stopCrashes(manifold::BCR4BP12Manifold)
    propagator = MBD.Propagator()
    crashEvent = DifferentialEquations.VectorContinuousCallback(primaryDistanceCondition, terminateAffect!, 3)
    manifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = Vector{MBD.BCR4BP12ManifoldArc}(undef, length(manifold.initialConditions))
    for a::Int64 = 1:length(manifold.initialConditions)
        arc::MBD.BCR4BP12Arc = propagateWithEvent(propagator, crashEvent, real(manifold.initialConditions[a]), [0.0, manifold.TOF], manifold.periodicOrbit.dynamicsModel)
        manifoldArcs[a] = MBD.BCR4BP12ManifoldArc(manifold.periodicOrbit, manifold.orbitTimes[a], manifold.ds[a], manifold.initialConditions[a], TOF = getTimeByIndex(arc, -1))
    end

    return manifoldArcs
end
