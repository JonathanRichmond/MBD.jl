"""
BCR4BP pseudo-manifold wrapper

Author: Jonathan Richmond
C: 6/17/25
"""

import DifferentialEquations
import MBD: BCR4BPPseudoManifold

export stopCrashes

"""
    stopCrashes(pseudoManifold)

Return pseudo-manifold arcs, stopping propagation when a primary is encountered

# Arguments
- `manifold::BCR4BPPseudoManifold`: BCR4BP pseudo-manifold object
"""
function stopCrashes(pseudoManifold::BCR4BPPseudoManifold)
    propagator = MBD.Propagator()
    crashEvent = DifferentialEquations.VectorContinuousCallback(primaryDistanceCondition3, terminateAffectIndex!, 3)
    EarthRadius::Float64 = pseudoManifold.dynamicsModel.systemData.primaryData[1].bodyRadius/get12CharLength(pseudoManifold.dynamicsModel)
    MoonRadius::Float64 = pseudoManifold.dynamicsModel.systemData.primaryData[2].bodyRadius/get12CharLength(pseudoManifold.dynamicsModel)
    SunRadius::Float64 = pseudoManifold.dynamicsModel.systemData.primaryData[3].bodyRadius/get12CharLength(pseudoManifold.dynamicsModel)
    pseudoManifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc} = Vector{MBD.BCR4BPPseudoManifoldArc}(undef, length(pseudoManifold.initialConditions))
    for a::Int64 = 1:length(pseudoManifold.initialConditions)
        arc::MBD.BCR4BP12Arc = propagateWithEvent(propagator, crashEvent, pseudoManifold.initialConditions[a], [0.0, pseudoManifold.TOF], pseudoManifold.dynamicsModel, [pseudoManifold.dynamicsModel, EarthRadius, MoonRadius, SunRadius])
        pseudoManifoldArcs[a] = MBD.BCR4BPPseudoManifoldArc(pseudoManifold.periodicOrbit, pseudoManifold.dynamicsModel, pseudoManifold.orbitTimes[a], pseudoManifold.theta40, pseudoManifold.initialConditions[a], getTimeByIndex(arc, -1))
    end

    return pseudoManifoldArcs
end
