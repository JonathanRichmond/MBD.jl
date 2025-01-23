"""
State match constraint wrapper

Author: Jonathan Richmond
C: 9/23/22
U: 1/16/25
"""

import MBD: StateMatchConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(stateMatchConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `stateMatchConstraint::StateMatchConstraint`: State match constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(stateMatchConstraint::StateMatchConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    return getData(stateMatchConstraint.variable1)[stateMatchConstraint.constrainedIndices]-getData(stateMatchConstraint.variable2)[stateMatchConstraint.constrainedIndices]
end

"""
    getNumConstraintRows(stateMatchConstraint)

Return number of constraints

# Arguments
- `stateMatchConstraint::StateMatchConstraint`: State match constraint object
"""
function getNumConstraintRows(stateMatchConstraint::StateMatchConstraint)
    return length(stateMatchConstraint.constrainedIndices)
end

"""
    getPartials_ConstraintWRTVariables(stateMatchConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `stateMatchConstraint::StateMatchConstraint`: State match constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(stateMatchConstraint::StateMatchConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    numConstraints::Int16 = Int16(getNumConstraintRows(stateMatchConstraint))
    partials1::Matrix{Float64} = zeros(Float64, (numConstraints,length(stateMatchConstraint.variable1.data)))
    partials2::Matrix{Float64} = zeros(Float64, (numConstraints,length(stateMatchConstraint.variable2.data)))
    [(partials1[r,stateMatchConstraint.constrainedIndices[r]] = 1) for r in Int16(1):numConstraints]
    [(partials2[r,stateMatchConstraint.constrainedIndices[r]] = -1) for r in Int16(1):numConstraints]
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}(stateMatchConstraint.variable1 => partials1, stateMatchConstraint.variable2 => partials2)

    return partialsMap
end

"""
    shallowClone(stateMatchConstraint)

Return copy of state match constraint object

# Arguments
- `stateMatchConstraint::StateMatchConstraint`: State match constraint object
"""
function shallowClone(stateMatchConstraint::StateMatchConstraint)
    return StateMatchConstraint(stateMatchConstraint.variable1, stateMatchConstraint.variable2, [Int64(i) for i in stateMatchConstraint.constrainedIndices])
end

"""
    updatePointers!(stateMatchConstraint, copiedObjectMap)

Update pointers for state match constraint object

# Arguments
- `stateMatchConstraint::StateMatchConstraint`: State match constraint object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(stateMatchConstraint::StateMatchConstraint, copiedObjectMap::Dict)
    stateMatchConstraint.variable1 = updatePointer(stateMatchConstraint.variable1, copiedObjectMap, true)
    stateMatchConstraint.variable2 = updatePointer(stateMatchConstraint.variable2, copiedObjectMap, true)
end
