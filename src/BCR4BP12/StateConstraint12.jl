"""
BCR4BP P1-P2 state constraint wrapper

Author: Jonathan Richmond
C: 4/9/25
U: 4/15/25
"""

import MBD: BCR4BP12StateConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(stateConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `stateConstraint::BCR4BP12StateConstraint`: BCR4BP P1-P2 state constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(stateConstraint::BCR4BP12StateConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    return getData(stateConstraint.variable)[stateConstraint.constrainedIndices]-stateConstraint.values
end

"""
    getNumConstraintRows(stateConstraint)

Return number of constraints

# Arguments
- `stateConstraint::BCR4BP12StateConstraint`: BCR4BP P1-P2 state constraint object
"""
function getNumConstraintRows(stateConstraint::BCR4BP12StateConstraint)
    return length(stateConstraint.constrainedIndices)
end

"""
    getPartials_ConstraintWRTVariables(stateConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `stateConstraint::BCR4BP12StateConstraint`: BCR4BP P1-P2 state constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(stateConstraint::BCR4BP12StateConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    numConstraints::Int16 = Int16(getNumConstraintRows(stateConstraint))
    partials::Matrix{Float64} = zeros(Float64, (numConstraints,length(stateConstraint.variable.data)))
    [(partials[r,stateConstraint.constrainedIndices[r]] = 1) for r in Int16(1):numConstraints]
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}(stateConstraint.variable => partials)

    return partialsMap
end

"""
    shallowClone(stateConstraint, dynamicsModel)

Return copy of state constraint object

# Arguments
- `stateConstraint::BCR4BP12StateConstraint`: BCR4BP P1-P2 state constraint object
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
function shallowClone(stateConstraint::BCR4BP12StateConstraint, dynamicsModel::BCR4BP12DynamicsModel)
    node = MBD.BCR4BP12Node(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    object = BCR4BP12StateConstraint(node, [Int64(i) for i in stateConstraint.constrainedIndices], stateConstraint.values)
    object.constrainedIndices = stateConstraint.constrainedIndices
    object.values = stateConstraint.values
    object.variable = stateConstraint.variable
    
    return object
end

"""
    updatePointers!(stateConstraint, copiedObjectMap)

Update pointers for state constraint object

# Arguments
- `stateConstraint::BCR4BP12StateConstraint`: BCR4BP P1-P2 state constraint object
- `copiedObjectMap::IdDict{Any, Any}`: Map between old and new objects
"""
function updatePointers!(stateConstraint::BCR4BP12StateConstraint, copiedObjectMap::IdDict{Any, Any})
    stateConstraint.variable = updatePointer(stateConstraint.variable, copiedObjectMap, true)
end
