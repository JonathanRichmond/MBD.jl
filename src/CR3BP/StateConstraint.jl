"""
CR3BP state constraint wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 1/16/25
"""

import MBD: CR3BPStateConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(stateConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `stateConstraint::CR3BPStateConstraint`: CR3BP state constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(stateConstraint::CR3BPStateConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    return getData(stateConstraint.variable)[stateConstraint.constrainedIndices]-stateConstraint.values
end

"""
    getNumConstraintRows(stateConstraint)

Return number of constraints

# Arguments
- `stateConstraint::CR3BPStateConstraint`: CR3BP state constraint object
"""
function getNumConstraintRows(stateConstraint::CR3BPStateConstraint)
    return length(stateConstraint.constrainedIndices)
end

"""
    getPartials_ConstraintWRTVariables(stateConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `stateConstraint::CR3BPStateConstraint`: CR3BP state constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(stateConstraint::CR3BPStateConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    numConstraints::Int16 = Int16(getNumConstraintRows(stateConstraint))
    partials::Matrix{Float64} = zeros(Float64, (numConstraints,length(stateConstraint.variable.data)))
    [(partials[r,stateConstraint.constrainedIndices[r]] = 1) for r in Int16(1):numConstraints]
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}(stateConstraint.variable => partials)

    return partialsMap
end

"""
    shallowClone(stateConstraint)

Return copy of state constraint object

# Arguments
- `stateConstraint::CR3BPStateConstraint`: CR3BP state constraint object
"""
function shallowClone(stateConstraint::CR3BPStateConstraint)
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    node = MBD.CR3BPNode(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    object = CR3BPStateConstraint(node, [Int64(i) for i in stateConstraint.constrainedIndices], stateConstraint.values)
    object.constrainedIndices = stateConstraint.constrainedIndices
    object.values = stateConstraint.values
    object.variable = stateConstraint.variable
    
    return object
end

"""
    updatePointers!(stateConstraint, copiedObjectMap)

Update pointers for state constraint object

# Arguments
- `stateConstraint::CR3BPStateConstraint`: CR3BP state constraint object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(stateConstraint::CR3BPStateConstraint, copiedObjectMap::Dict)
    stateConstraint.variable = updatePointer(stateConstraint.variable, copiedObjectMap, true)
end
