"""
Segment wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 8/6/23
"""

import MBD: Segment

export getFinalState!, getFinalStateRate!, getPartials_FinalStateWRTEpoch!
export getPartials_FinalStateWRTInitialState!, getPartials_FinalStateWRTParams!
export getPropagatorParametersData, getVariables, lazyPropagate!, propagate!
export resetPropagatedArc

"""
    getFinalState!(segment)

Return final state

# Arguments
- `segment::Segment`: Segment object
"""
function getFinalState!(segment::Segment)
    lazyPropagate!(segment, MBD.SIMPLE)

    return getStateByIndex(segment.propArc, -1)
end

"""
    getFinalStateRate!(segment)

Return time derivative of final state

# Arguments
- `segment::Segment`: Segment object
"""
function getFinalStateRate!(segment::Segment)
    lazyPropagate!(segment, MBD.SIMPLE)

    return evaluateEquations(segment.propArc.dynamicsModel, MBD.SIMPLE, getTimeByIndex(segment.propArc, -1), getStateByIndex(segment.propArc, -1), getPropagatorParametersData(segment))
end

"""
    getPartials_FinalStateWRTEpoch!(segment)

Return partial derivatives of final state with respect to origin epoch

# Arguments
- `segment::Segment`: Segment object
"""
function getPartials_FinalStateWRTEpoch!(segment::Segment)
    lazyPropagate!(segment, MBD.FULL)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)
    dqdT::Vector{Float64} = getEpochDependencies(segment.propArc.dynamicsModel, qf)

    return reshape(dqdT, :, 1)
end

"""
    getPartials_FinalStateWRTInitialState(segment)

Return partial derivatives of final state with respect to initial state

# Arguments
- `segment::Segment`: Segment object
"""
function getPartials_FinalStateWRTInitialState!(segment::Segment)
    lazyPropagate!(segment, MBD.STM)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)

    return getStateTransitionMatrix(segment.propArc.dynamicsModel, qf)
end

"""
    getPartials_FinalStateWRTParams!(segment)

Return partial derivatives of final state with respect to parameters

# Arguments
- `segment::Segment`: Segment object
"""
function getPartials_FinalStateWRTParams!(segment::Segment)
    lazyPropagate!(segment, MBD.FULL)
    qf::Vector{Float64} = getStateByIndex(segment.propArc, -1)

    return getParameterDependencies(segment.propArc.dynamicsModel, qf)
end

"""
    getPropagatorParametersData(segment)

Return propagation parameters

# Arguments
- `segment::Segment`: Segment object
"""
function getPropagatorParametersData(segment::Segment)
    return getData(segment.propParams)
end

"""
    getVariables(segment)

Return variables

# Arguments
- `segment::Segment`: Segment object
"""
function getVariables(segment::Segment)
    return [segment.TOF, segment.propParams]
end

"""
    lazyPropagate!(segment, minEquationType)

Return segment object propagated with minimum required EOMs

# Arguments
- `segment::Segment`: Segment object
- `minEquationType::EquationType`: EOM type
"""
function lazyPropagate!(segment::Segment, minEquationType::MBD.EquationType)
    isempty(segment.propArc.times) && propagate!(segment, minEquationType)
    q::Vector{Float64} = getStateByIndex(segment.propArc, 1)
    (length(q) < getStateSize(segment.propArc.dynamicsModel, minEquationType)) && propagate!(segment, minEquationType)
end

"""
    propagate!(segment, equationType)

Return propagated segment

# Arguments
- `segment::Segment`: Segment object
- `equationType::EquationType`: EOM type
"""
function propagate!(segment::Segment, equationType::MBD.EquationType)
    q0::Vector{Float64} = appendExtraInitialConditions(segment.originNode.dynamicsModel, getData(segment.originNode.state), equationType)
    t0::Float64 = getData(segment.originNode.epoch)[1]
    tSpan::Vector{Float64} = [t0, t0+getData(segment.TOF)[1]]
    segment.propagator.equationType = equationType
    segment.propArc = propagate(segment.propagator, q0, tSpan, segment.originNode.dynamicsModel, getPropagatorParametersData(segment))
end

"""
    resetPropagatedArc!(segment)

Return empty arc object

# Arguments
- `segment::Segment`: Segment object
"""
function resetPropagatedArc!(segment::Segment)
    segment.propArc = MBD.Arc(segment.originNode.dynamicsModel)
end

"""
    shallowClone(segment)

Return copy of segment object

# Arguments
- `segment::Segment`: Segment object
"""
function shallowClone(segment::Segment)
    object = Segment(segment.TOF.data[1], segment.originNode, segment.terminalNode)
    object.originNode = segment.originNode
    object.propArc = segment.propArc
    object.propagator = segment.propagator
    object.propParams = segment.propParams
    object.terminalNode = segment.terminalNode
    object.TOF = segment.TOF

    return object
end

"""
    updatePointers!(segment, copiedObjectMap)

Update pointers for segment object

# Arguments
- `segment::Segment`: Segment object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(segment::Segment, copiedObjectMap::Dict)
    segment.TOF = updatePointer(segment.TOF, copiedObjectMap, true)
    segment.propParams = updatePointer(segment.propParams, copiedObjectMap, true)
    segment.originNode = updatePointer(segment.originNode, copiedObjectMap, true)
    segment.terminalNode = updatePointer(segment.terminalNode, copiedObjectMap, true)
end
