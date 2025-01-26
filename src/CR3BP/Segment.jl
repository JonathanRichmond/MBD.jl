"""
Segment wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 1/15/25
"""

import MBD: CR3BPSegment

export getFinalState!, getFinalStateRate!, getPartials_FinalStateWRTEpoch!
export getPartials_FinalStateWRTInitialState!, getVariables, lazyPropagate!, propagate!
export resetPropagatedArc!

"""
    getFinalState!(segment)

Return final state

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
function getFinalState!(segment::CR3BPSegment)
    lazyPropagate!(segment, MBD.SIMPLE)

    return getStateByIndex(segment.propArc, -1)
end

"""
    getFinalStateRate!(segment)

Return time derivative of final state

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
function getFinalStateRate!(segment::CR3BPSegment)
    lazyPropagate!(segment, MBD.SIMPLE)

    return evaluateEquations(segment.propArc.dynamicsModel, MBD.SIMPLE, getTimeByIndex(segment.propArc, -1), getStateByIndex(segment.propArc, -1))
end

"""
    getPartials_FinalStateWRTEpoch!(segment)

Return partial derivatives of final state with respect to origin epoch

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
function getPartials_FinalStateWRTEpoch!(segment::CR3BPSegment)
    lazyPropagate!(segment, MBD.FULL)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)
    dqdT::Vector{Float64} = getEpochDependencies(segment.propArc.dynamicsModel, qf)

    return reshape(dqdT, :, 1)
end

"""
    getPartials_FinalStateWRTInitialState(segment)

Return partial derivatives of final state with respect to initial state

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
function getPartials_FinalStateWRTInitialState!(segment::CR3BPSegment)
    lazyPropagate!(segment, MBD.STM)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)

    return getStateTransitionMatrix(segment.propArc.dynamicsModel, qf)
end

"""
    getVariables(segment)

Return variables

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
function getVariables(segment::CR3BPSegment)
    return [segment.TOF]
end

"""
    lazyPropagate!(segment, minEquationType)

Return segment object propagated with minimum required EOMs

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
- `minEquationType::EquationType`: EOM type
"""
function lazyPropagate!(segment::CR3BPSegment, minEquationType::MBD.EquationType)
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
- `segment::CR3BPSegment`: CR3BP segment object
- `equationType::EquationType`: EOM type
"""
function propagate!(segment::CR3BPSegment, equationType::MBD.EquationType)
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
- `segment::CR3BPSegment`: CR3BP segment object
"""
function resetPropagatedArc!(segment::CR3BPSegment)
    segment.propArc = MBD.CR3BPArc(segment.originNode.dynamicsModel)
end

"""
    shallowClone(segment)

Return copy of segment object

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
function shallowClone(segment::CR3BPSegment)
    object = CR3BPSegment(segment.TOF.data[1], segment.originNode, segment.terminalNode)
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
- `segment::CR3BPSegment`: CR3BP segment object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(segment::CR3BPSegment, copiedObjectMap::Dict)
    segment.TOF = updatePointer(segment.TOF, copiedObjectMap, true)
    segment.originNode = updatePointer(segment.originNode, copiedObjectMap, true)
    segment.terminalNode = updatePointer(segment.terminalNode, copiedObjectMap, true)
end
