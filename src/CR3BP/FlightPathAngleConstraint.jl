"""
CR3BP flight path angle constraint wrapper

Author: Jonathan Richmond
C: 7/1/25
U: 7/2/25
"""

import MBD: CR3BPFlightPathAngleConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(flightPathAngleConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `flightPathAngleConstraint::CR3BPFlightPathAngleConstraint`: CR3BP flight path angle constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(flightPathAngleConstraint::CR3BPFlightPathAngleConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    primaryState::Vector{Float64} = getPrimaryState(flightPathAngleConstraint.dynamicsModel, flightPathAngleConstraint.primary)
    q_PC::Vector{Float64} = flightPathAngleConstraint.variable.data[1:6]-primaryState
    
    return [LinearAlgebra.dot(q_PC[1:3], q_PC[4:6])/(LinearAlgebra.norm(q_PC[1:3])*LinearAlgebra.norm(q_PC[4:6]))-sin(flightPathAngleConstraint.value)]
end

"""
    getNumConstraintRows(flightPathAngleConstraint)

Return number of constraints

# Arguments
- `flightPathAngleConstraint::CR3BPFlightPathAngleConstraint`: CR3BP flight path angle constraint object
"""
function getNumConstraintRows(flightPathAngleConstraint::CR3BPFlightPathAngleConstraint)
    return 1
end

"""
    getPartials_ConstraintWRTVariables(flightPathAngleConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `flightPathAngleConstraint::CR3BPFlightPathAngleConstraint`: CR3BP flight path angle constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(flightPathAngleConstraint::CR3BPFlightPathAngleConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    numConstraints::Int16 = Int16(getNumConstraintRows(flightPathAngleConstraint))
    partials::Matrix{Float64} = zeros(Float64, (numConstraints,length(flightPathAngleConstraint.variable.data)))

    primaryState::Vector{Float64} = getPrimaryState(flightPathAngleConstraint.dynamicsModel, flightPathAngleConstraint.primary)
    q_PC::Vector{Float64} = flightPathAngleConstraint.variable.data[1:6]-primaryState
    [partials[1,j] = (q_PC[j+3]*LinearAlgebra.norm(q_PC[1:3])^(2)-LinearAlgebra.dot(q_PC[1:3], q_PC[4:6])*q_PC[j])/(LinearAlgebra.norm(q_PC[1:3])^(3)*LinearAlgebra.norm(q_PC[4:6])) for j = 1:3]
    [partials[1,j] = (q_PC[j-3]*LinearAlgebra.norm(q_PC[4:6])^(2)-LinearAlgebra.dot(q_PC[1:3], q_PC[4:6])*q_PC[j])/(LinearAlgebra.norm(q_PC[1:3])*LinearAlgebra.norm(q_PC[4:6])^(3)) for j = 4:6]
    
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}(flightPathAngleConstraint.variable => partials)

    return partialsMap
end

"""
    shallowClone(flightPathAngleConstraint, dynamicsModel)

Return copy of flight path angle constraint object

# Arguments
- `flightPathAngleConstraint::CR3BPFlightPathAngleConstraint`: CR3BP flight path angle constraint object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
function shallowClone(flightPathAngleConstraint::CR3BPFlightPathAngleConstraint, dynamicsModel::MBD.CR3BPDynamicsModel)
    node = MBD.CR3BPNode(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    object = CR3BPFlightPathAngleConstraint(node, flightPathAngleConstraint.primary, flightPathAngleConstraint.value)
    object.dynamicsModel = flightPathAngleConstraint.dynamicsModel
    object.primary = flightPathAngleConstraint.primary
    object.value = flightPathAngleConstraint.value
    object.variable = flightPathAngleConstraint.variable
    
    return object
end

"""
    updatePointers!(flightPathAngleConstraint, copiedObjectMap)

Update pointers for flight path angle constraint object

# Arguments
- `flightPathAngleConstraint::CR3BPFlightPathAngleConstraint`: CR3BP flight path angle constraint object
- `copiedObjectMap::IdDict{Any, Any}`: Map between old and new objects
"""
function updatePointers!(flightPathAngleConstraint::CR3BPFlightPathAngleConstraint, copiedObjectMap::IdDict{Any, Any})
    flightPathAngleConstraint.variable = updatePointer(flightPathAngleConstraint.variable, copiedObjectMap, true)
end
