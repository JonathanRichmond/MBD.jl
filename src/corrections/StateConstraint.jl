"""
State constraint wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 8/22/23
"""

import MBD: StateConstraint

export evaluateConstraint, getNumberConstraintRows
export getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(stateConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `stateConstraint::StateConstraint`: State constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(stateConstraint::StateConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int64}, freeVariableVector::Vector{Float64})
    data::Vector{Float64} = getData(stateConstraint.variable)

    return data[stateConstraint.constrainedIndices]-stateConstraint.values
end

"""
    getNumberConstraintRows(stateConstraint)

Return number of constraints

# Arguments
- `stateConstraint::StateConstraint`: State constraint object
"""
function getNumberConstraintRows(stateConstraint::StateConstraint)
    return length(stateConstraint.constrainedIndices)
end

"""
    getPartials_ConstraintWRTVariables(stateConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `stateConstraint::StateConstraint`: State constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(stateConstraint::StateConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int64}, freeVariableVector::Vector{Float64})
    partials::Matrix{Float64} = zeros(Float64, (length(stateConstraint.constrainedIndices), length(stateConstraint.variable.data)))
    [(partials[r,stateConstraint.constrainedIndices[r]] = 1) for r in 1:length(stateConstraint.constrainedIndices)]
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}()
    partialsMap[stateConstraint.variable] = partials

    return partialsMap
end

"""
    shallowClone(stateConstraint)

Return copy of state constraint object

# Arguments
- `stateConstraint::StateConstraint`: State constraint object
"""
function shallowClone(stateConstraint::StateConstraint)
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    node = Node(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    object = StateConstraint(node, stateConstraint.constrainedIndices, stateConstraint.values)
    object.constrainedIndices = stateConstraint.constrainedIndices
    object.values = stateConstraint.values
    object.variable = stateConstraint.variable
    
    return object
end

"""
    updatePointers!(stateConstraint, copiedObjectMap)

Update pointers for state constraint object

# Arguments
- `stateConstraint::StateConstraint`: State constraint object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(stateConstraint::StateConstraint, copiedObjectMap::Dict)
    stateConstraint.variable = updatePointer(stateConstraint.variable, copiedObjectMap, true)
end
