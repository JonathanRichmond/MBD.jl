"""
BCR4BP P1-P2 continuity constraint wrapper

Author: Jonathan Richmond
C: 4/9/25
"""

import StaticArrays
import MBD: BCR4BP12ContinuityConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(continuityConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `continuityConstraint::BCR4BP12ContinuityConstraint`: BCR4BP P1-P2 continuity constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(continuityConstraint::BCR4BP12ContinuityConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    propState::Vector{Float64} = getFinalState!(continuityConstraint.segment)
    terminalNodeState::StaticArrays.SVector{7, Float64} = StaticArrays.SVector{7, Float64}(getData(continuityConstraint.segment.terminalNode.state))
    
    return [propState[continuityConstraint.constrainedIndices[index]]-terminalNodeState[continuityConstraint.constrainedIndices[index]] for index in 1:getNumConstraintRows(continuityConstraint)]
end

"""
    getNumConstraintRows(continuityConstraint)

Return number of constraints

# Arguments
- `continuityConstraint::BCR4BP12ContinuityConstraint`: BCR4BP P1-P2 continuity constraint object
"""
function getNumConstraintRows(continuityConstraint::BCR4BP12ContinuityConstraint)
    return length(continuityConstraint.constrainedIndices)
end

"""
    getPartials_ConstraintWRTVariables(continuityConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `continuityConstraint::BCR4BP12ContinuityConstraint`: BCR4BP P1-P2 continuity constraint object
- `freeVariableIndexMap::Dict{Variable, Int16}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(continuityConstraint::BCR4BP12ContinuityConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    numConstraints::Int16 = Int16(getNumConstraintRows(continuityConstraint))
    STM::StaticArrays.SMatrix{7, 7, Float64} = StaticArrays.SMatrix{7, 7, Float64}(getPartials_FinalStateWRTInitialState!(continuityConstraint.segment))
    propStateRate::StaticArrays.SVector{7, Float64} = StaticArrays.SVector{7, Float64}(getFinalStateRate!(continuityConstraint.segment))
    finalStateWRTTime::Matrix{Float64} = zeros(Float64, (numConstraints,1))
    finalStateWRTInitialState::Matrix{Float64} = zeros(Float64, (numConstraints,6))
    finalStateWRTTargetState::Matrix{Float64} = copy(finalStateWRTInitialState)
    for index::Int16 in Int16(1):numConstraints
        finalStateWRTTime[index,1] = propStateRate[continuityConstraint.constrainedIndices[index]]
        finalStateWRTInitialState[index,:] = STM[continuityConstraint.constrainedIndices[index],:]
        finalStateWRTTargetState[index, continuityConstraint.constrainedIndices[index]] = -1
    end
    partials::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}()
    partials[continuityConstraint.segment.TOF] = finalStateWRTTime
    partials[continuityConstraint.segment.originNode.state] = finalStateWRTInitialState
    partials[continuityConstraint.segment.terminalNode.state] = finalStateWRTTargetState
    if !isEpochIndependent(continuityConstraint.segment.originNode.dynamicsModel)
        dqdT::StaticArrays.SMatrix{7, 1, Float64} = StaticArrays.SMatrix{7, 1, Float64}(getPartials_FinalStateWRTEpoch!(continuityConstraint.segment))
        finalStateWRTEpoch::Matrix{float64} = zeros(Float64, (numConstraints,1))
        [finalStateWRTEpoch[index,1] = dqdT[continuityConstraint.constrainedIndices[index],1] for index in Int16(1):numConstraints]
        partials[continuityConstraint.segment.originNode.epoch] = finalStateWRTEpoch
    end

    return partials
end

"""
    shallowClone(continuityConstraint, dynamicsModel)

Return copy of continuity constraint object

# Arguments
- `continuityConstraint::BCR4BP12ContinuityConstraint`: BCR4BP P1-P2 continuity constraint object
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
function shallowClone(continuityConstraint::BCR4BP12ContinuityConstraint, dynamicsModel::BCR4BP12DynamicsModel)
    object = BCR4BP12ContinuityConstraint(continuityConstraint.segment)
    object.constrainedIndices = continuityConstraint.constrainedIndices
    object.segment = continuityConstraint.segment
    
    return object
end

"""
    updatePointers!(continuityConstraint, copiedObjectMap)

Update pointers for continuity constraint object

# Arguments
- `continuityConstraint::BCR4BP12ContinuityConstraint`: BCR4BP P1-P2 continuity constraint object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(continuityConstraint::BCR4BP12ContinuityConstraint, copiedObjectMap::Dict)
    continuityConstraint.segment = updatePointer(continuityConstraint.segment, copiedObjectMap, true)
end
