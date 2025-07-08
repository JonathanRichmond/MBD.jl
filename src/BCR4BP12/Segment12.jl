"""
Segment wrapper

Author: Jonathan Richmond
C: 4/9/25
U: 7/8/25
"""

import MBD: BCR4BP12Segment

export getFinalState!, getFinalStateRate!, getPartials_FinalStateWRTEpoch!
export getPartials_FinalStateWRTInitialState!, getVariables, lazyPropagate!, propagate!
export resetPropagatedArc!, updateTerminalNodeEpoch!

"""
    getFinalState!(segment)

Return final state

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function getFinalState!(segment::BCR4BP12Segment)
    lazyPropagate!(segment, MBD.SIMPLE)

    return getStateByIndex(segment.propArc, -1)
end

"""
    getFinalStateRate!(segment)

Return time derivative of final state

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function getFinalStateRate!(segment::BCR4BP12Segment)
    lazyPropagate!(segment, MBD.SIMPLE)

    return evaluateEquations(segment.propArc.dynamicsModel, MBD.SIMPLE, getTimeByIndex(segment.propArc, -1), getStateByIndex(segment.propArc, -1))
end

"""
    getPartials_FinalStateWRTEpoch!(segment)

Return partial derivatives of final state with respect to origin epoch

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function getPartials_FinalStateWRTEpoch!(segment::BCR4BP12Segment)
    lazyPropagate!(segment, MBD.FULL)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)
    dqdT::Vector{Float64} = getEpochDependencies(segment.propArc.dynamicsModel, qf)

    return reshape(dqdT, :, 1)
end

"""
    getPartials_FinalStateWRTInitialState(segment)

Return partial derivatives of final state with respect to initial state

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function getPartials_FinalStateWRTInitialState!(segment::BCR4BP12Segment)
    lazyPropagate!(segment, MBD.STM)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)

    return getStateTransitionMatrix(segment.propArc.dynamicsModel, qf)
end

"""
    getVariables(segment)

Return variables

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function getVariables(segment::BCR4BP12Segment)
    return [segment.TOF]
end

"""
    lazyPropagate!(segment, minEquationType)

Return segment object propagated with minimum required EOMs

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
- `minEquationType::EquationType`: EOM type
"""
function lazyPropagate!(segment::BCR4BP12Segment, minEquationType::MBD.EquationType)
    if isempty(segment.propArc.times)
        propagate!(segment, minEquationType)
    elseif length(getStateByIndex(segment.propArc, 1)) < getStateSize(segment.propArc.dynamicsModel, minEquationType)
        propagate!(segment, minEquationType)
    end
end

"""
    propagate!(segment, equationType)

Return propagated segment

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
- `equationType::EquationType`: EOM type
"""
function propagate!(segment::BCR4BP12Segment, equationType::MBD.EquationType)
    q0::Vector{Float64} = appendExtraInitialConditions(segment.originNode.dynamicsModel, getData(segment.originNode.state), equationType)
    t0::Float64 = getData(segment.originNode.epoch)[1]
    tSpan::Vector{Float64} = [t0, t0+getData(segment.TOF)[1]]
    segment.propagator.equationType = equationType
    segment.propArc = propagate(segment.propagator, q0, tSpan, segment.originNode.dynamicsModel)
end

"""
    resetPropagatedArc!(segment)

Return empty arc object

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function resetPropagatedArc!(segment::BCR4BP12Segment)
    segment.propArc = MBD.BCR4BP12Arc(segment.originNode.dynamicsModel)
end

"""
    shallowClone(segment)

Return copy of segment object

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function shallowClone(segment::BCR4BP12Segment)
    object = BCR4BP12Segment(segment.TOF.data[1], segment.originNode, segment.terminalNode)
    object.originNode = segment.originNode
    object.propArc = segment.propArc
    object.propagator = segment.propagator
    object.terminalNode = segment.terminalNode
    object.TOF = segment.TOF

    return object
end

"""
    updatePointers!(segment, copiedObjectMap)

Update pointers for segment object

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
- `copiedObjectMap::IdDict{Any, Any}`: Map between old and new objects
"""
function updatePointers!(segment::BCR4BP12Segment, copiedObjectMap::IdDict{Any, Any})
    segment.TOF = updatePointer(segment.TOF, copiedObjectMap, true)
    segment.originNode = updatePointer(segment.originNode, copiedObjectMap, true)
    segment.terminalNode = updatePointer(segment.terminalNode, copiedObjectMap, true)
end

"""
    updateTerminalNodeEpoch!(segment)

Update epoch of terminal node

# Arguments
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function updateTerminalNodeEpoch!(segment::BCR4BP12Segment)
    segment.terminalNode.epoch.data = getData(segment.originNode.epoch)+getData(segment.TOF)
end
