"""
Continuity constraint wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 8/11/23
"""

import MBD: ContinuityConstraint

export evaluateConstraint, getNumberConstraintRows
export getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(continuityConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `continuityConstraint::ContinuityConstraint`: Continuity constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(continuityConstraint::ContinuityConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int64}, freeVariableVector::Vector{Float64})
    propState::Vector{Float64} = getFinalState!(continuityConstraint.segment)
    terminalNodeState::Vector{Float64} = getData(continuityConstraint.segment.terminalNode.state)
    
    return [propState[continuityConstraint.constrainedIndices[index]]-terminalNodeState[continuityConstraint.constrainedIndices[index]] for index in 1:length(continuityConstraint.constrainedIndices)]
end

"""
    getNumberConstraintRows(continuityConstraint)

Return number of constraints

# Arguments
- `continuityConstraint::ContinuityConstraint`: Continuity constraint object
"""
function getNumberConstraintRows(continuityConstraint::ContinuityConstraint)
    return length(continuityConstraint.constrainedIndices)
end

"""
    getPartials_ConstraintWRTVariables(continuityConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `continuityConstraint::ContinuityConstraint`: Continuity constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(continuityConstraint::ContinuityConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int64}, freeVariableVector::Vector{Float64})
    STM::Matrix{Float64} = getPartials_FinalStateWRTInitialState!(continuityConstraint.segment)
    propStateRate::Vector{Float64} = getFinalStateRate!(continuityConstraint.segment)
    finalStateWRTTime::Matrix{Float64} = zeros(Float64, (length(continuityConstraint.constrainedIndices), 1))
    finalStateWRTInitialState::Matrix{Float64} = zeros(Float64, (length(continuityConstraint.constrainedIndices), size(STM, 2)))
    finalStateWRTTargetState::Matrix{Float64} = copy(finalStateWRTInitialState)
    for index::Int64 in 1:length(continuityConstraint.constrainedIndices)
        finalStateWRTTime[index,1] = propStateRate[continuityConstraint.constrainedIndices[index]]
        finalStateWRTInitialState[index,:] = STM[continuityConstraint.constrainedIndices[index],:]
        finalStateWRTTargetState[index, continuityConstraint.constrainedIndices[index]] = -1
    end
    partials::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}()
    partials[continuityConstraint.segment.TOF] = finalStateWRTTime
    partials[continuityConstraint.segment.originNode.state] = finalStateWRTInitialState
    partials[continuityConstraint.segment.terminalNode.state] = finalStateWRTTargetState
    if !isEpochIndependent(continuityConstraint.segment.originNode.dynamicsModel)
        dqdT::Matrix{Float64} = getPartials_FinalStateWRTEpoch!(continuityConstraint.segment)
        finalStateWRTEpoch::Matrix{float64} = zeros(Float64, (length(continuityConstraint.constrainedIndices), 1))
        [finalStateWRTEpoch[index,1] = dqdT[continuityConstraint.constrainedIndices[index],1] for index in 1:length(continuityConstraint.constrainedIndices)]
        partials[continuityConstraint.segment.originNode.epoch] = finalStateWRTEpoch
    end
    if !isempty(continuityConstraint.segment.propParams.data)
        dqdp::Matrix{Float64} = getPartials_FinalStateWRTParams!(continuityConstraint.segment)
        finalStateWRTParams::Matrix{Float64} = zeros(Float64, (length(continuityConstraint.constrainedIndices), size(dqdp, 2)))
        [finalStateWRTParams[index] = dqdp[continuityConstraint.constrainedIndices[index]] for index in 1:length(continuityConstraint.constrainedIndices)]
        partials[getPropagatorParametersData(continuityConstraint.segment)] = finalStateWRTParams
    end

    return partials
end

"""
    shallowClone(continuityConstraint)

Return copy of continuity constraint object

# Arguments
- `continuityConstraint::ContinuityConstraint`: Continuity constraint object
"""
function shallowClone(continuityConstraint::ContinuityConstraint)
    object = ContinuityConstraint(continuityConstraint.segment)
    object.constrainedIndices = continuityConstraint.constrainedIndices
    object.segment = continuityConstraint.segment
    
    return object
end
