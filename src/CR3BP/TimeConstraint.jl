"""
CR3BP time constraint wrapper

Author: Jonathan Richmond
C: 6/30/25
"""

import MBD: CR3BPTimeConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(timeConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `timeConstraint::CR3BPTimeConstraint`: CR3BP time constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(timeConstraint::CR3BPTimeConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    return getData(timeConstraint.variable)[1]-timeConstraint.value
end

"""
    getNumConstraintRows(timeConstraint)

Return number of constraints

# Arguments
- `timeConstraint::CR3BPTimeConstraint`: CR3BP time constraint object
"""
function getNumConstraintRows(timeConstraint::CR3BPTimeConstraint)
    return 1
end

"""
    getPartials_ConstraintWRTVariables(timeConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `timeConstraint::CR3BPTimeConstraint`: CR3BP time constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(timeConstraint::CR3BPTimeConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    partials::Matrix{Float64} = reshape([1], (1,1))
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}(timeConstraint.variable => partials)

    return partialsMap
end

"""
    shallowClone(timeConstraint, dynamicsModel)

Return copy of time constraint object

# Arguments
- `timeConstraint::CR3BPTimeConstraint`: CR3BP time constraint object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
function shallowClone(timeConstraint::CR3BPTimeConstraint, dynamicsModel::MBD.CR3BPDynamicsModel)
    node1 = MBD.CR3BPNode(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    node2 = MBD.CR3BPNode(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    segment = MBD.CR3BPSegment(0.0, node1, node2)
    object = CR3BPTimeConstraint(segment, timeConstraint.value)
    object.value = timeConstraint.value
    object.variable = timeConstraint.variable
    
    return object
end

"""
    updatePointers!(timeConstraint, copiedObjectMap)

Update pointers for time constraint object

# Arguments
- `timeConstraint::CR3BPTimeConstraint`: CR3BP time constraint object
- `copiedObjectMap::IdDict{Any, Any}`: Map between old and new objects
"""
function updatePointers!(timeConstraint::CR3BPTimeConstraint, copiedObjectMap::IdDict{Any, Any})
    timeConstraint.variable = updatePointer(timeConstraint.variable, copiedObjectMap, true)
end
