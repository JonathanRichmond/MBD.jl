"""
Jacobi constraint wrapper

Author: Jonathan Richmond
C: 9/23/22
U: 9/6/23
"""

import MBD: JacobiConstraint

export evaluateConstraint, getNumberConstraintRows
export getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(jacobiConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `jacobiConstraint::JacobiConstraint`: Jacobi constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(jacobiConstraint::JacobiConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int64}, freeVariableVector::Vector{Float64})
    return [getJacobiConstant(jacobiConstraint.dynamicsModel, getData(jacobiConstraint.state))-jacobiConstraint.value]
end

"""
    getNumberConstraintRows(jacobiConstraint)

Return number of constraints

# Arguments
- `jacobiConstraint::JacobiConstraint`: Jacobi constraint object
"""
function getNumberConstraintRows(jacobiConstraint::JacobiConstraint)
    return 1
end

"""
    getPartials_ConstraintWRTVariables(jacobiConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `jacobiConstraint::JacobiConstraint`: Jacobi constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(jacobiConstraint::JacobiConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int64}, freeVariableVector::Vector{Float64})
    mu::Float64 = getMassRatio(jacobiConstraint.dynamicsModel.systemData)
    q::Vector{Float64} = getData(jacobiConstraint.state)
    r_13::Float64 = sqrt((q[1]+mu)^2+q[2]^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-1+mu)^2+q[2]^2+q[3]^2)
    r_13_3::Float64 = r_13^3
    r_23_3::Float64 = r_23^3
    partialWRTNodeState::Matrix{Float64} = zeros(Float64, (1,6))
    partialWRTNodeState[1,:] = [-2*(q[1]+mu)*(1-mu)/r_13_3-2*mu*(q[1]-1+mu)/r_23_3+2*q[1], -2*q[2]*(1-mu)/r_13_3-2*q[2]*mu/r_23_3+2*q[2], -2*q[3]*(1-mu)/r_13_3-2*q[3]*mu/r_23_3, -2*q[4], -2*q[5], -2*q[6]]
    partials::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}()
    partials[jacobiConstraint.state] = partialWRTNodeState

    return partials
end

"""
    shallowClone(jacobiConstraint)

Return copy of Jacobi constraint object

# Arguments
- `jacobiConstraint::JacobiConstraint`: Jacobi constraint object
"""
function shallowClone(jacobiConstraint::JacobiConstraint)
    node = Node(jacobiConstraint.epoch.data[1], jacobiConstraint.state.data, jacobiConstraint.dynamicsModel)
    object = JacobiConstraint(node, jacobiConstraint.value)
    object.epoch = jacobiConstraint.epoch
    object.state = jacobiConstraint.state
    object.dynamicsModel = jacobiConstraint.dynamicsModel
    
    return object
end

"""
    updatePointers!(jacobiConstraint, copiedObjectMap)

Update pointers for Jacobi constraint object

# Arguments
- `jacobiConstraint::JacobiConstraint`: Jacobi constraint object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(jacobiConstraint::JacobiConstraint, copiedObjectMap::Dict)
    jacobiConstraint.epoch = updatePointer(jacobiConstraint.epoch, copiedObjectMap, true)
    jacobiConstraint.state = updatePointer(jacobiConstraint.state, copiedObjectMap, true)
end
