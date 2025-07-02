"""
CR3BP altitude constraint wrapper

Author: Jonathan Richmond
C: 7/1/25
U: 7/2/25
"""

import MBD: CR3BPAltitudeConstraint

export evaluateConstraint, getNumConstraintRows, getPartials_ConstraintWRTVariables

"""
    evaluateConstraint(altitudeConstraint, freeVariableIndexMap, freeVariableVector)

Return constraint error

# Arguments
- `altitudeConstraint::CR3BPAltitudeConstraint`: CR3BP altitude constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function evaluateConstraint(altitudeConstraint::CR3BPAltitudeConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    primaryPos::Vector{Float64} = getPrimaryState(altitudeConstraint.dynamicsModel, altitudeConstraint.primary)[1:3]
    primaryRad::Float64 = altitudeConstraint.dynamicsModel.systemData.primaryData[altitudeConstraint.primary].bodyRadius/getCharLength(altitudeConstraint.dynamicsModel)

    return [LinearAlgebra.norm(altitudeConstraint.variable.data[1:3]-primaryPos)-primaryRad-altitudeConstraint.value]
end

"""
    getNumConstraintRows(altitudeConstraint)

Return number of constraints

# Arguments
- `altitudeConstraint::CR3BPAltitudeConstraint`: CR3BP altitude constraint object
"""
function getNumConstraintRows(altitudeConstraint::CR3BPAltitudeConstraint)
    return 1
end

"""
    getPartials_ConstraintWRTVariables(altitudeConstraint, freeVariableIndexMap, freeVariableVector)

Return partial derivatives of constraint with respect to free variables

# Arguments
- `altitudeConstraint::CR3BPAltitudeConstraint`: CR3BP altitude constraint object
- `freeVariableIndexMap::Dict{Variable, Int64}`: Free variable index map
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function getPartials_ConstraintWRTVariables(altitudeConstraint::CR3BPAltitudeConstraint, freeVariableIndexMap::Dict{MBD.Variable, Int16}, freeVariableVector::Vector{Float64})
    numConstraints::Int16 = Int16(getNumConstraintRows(altitudeConstraint))
    partials::Matrix{Float64} = zeros(Float64, (numConstraints,length(altitudeConstraint.variable.data)))

    primaryPos::Vector{Float64} = getPrimaryState(altitudeConstraint.dynamicsModel, altitudeConstraint.primary)[1:3]
    d::Float64 = LinearAlgebra.norm(altitudeConstraint.variable.data[1:3]-primaryPos)
    [partials[1,j] = (altitudeConstraint.variable.data[j]-primaryPos[j])/d for j = 1:3]
    
    partialsMap::Dict{MBD.Variable, Matrix{Float64}} = Dict{MBD.Variable, Matrix{Float64}}(altitudeConstraint.variable => partials)

    return partialsMap
end

"""
    shallowClone(altitudeConstraint, dynamicsModel)

Return copy of altitude constraint object

# Arguments
- `altitudeConstraint::CR3BPAltitudeConstraint`: CR3BP altitude constraint object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
function shallowClone(altitudeConstraint::CR3BPAltitudeConstraint, dynamicsModel::MBD.CR3BPDynamicsModel)
    node = MBD.CR3BPNode(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dynamicsModel)
    object = CR3BPAltitudeConstraint(node, altitudeConstraint.primary, altitudeConstraint.value)
    object.dynamicsModel = altitudeConstraint.dynamicsModel
    object.primary = altitudeConstraint.primary
    object.value = altitudeConstraint.value
    object.variable = altitudeConstraint.variable
    
    return object
end

"""
    updatePointers!(altitudeConstraint, copiedObjectMap)

Update pointers for altitude constraint object

# Arguments
- `altitudeConstraint::CR3BPAltitudeConstraint`: CR3BP altitude constraint object
- `copiedObjectMap::IdDict{Any, Any}`: Map between old and new objects
"""
function updatePointers!(altitudeConstraint::CR3BPAltitudeConstraint, copiedObjectMap::IdDict{Any, Any})
    altitudeConstraint.variable = updatePointer(altitudeConstraint.variable, copiedObjectMap, true)
end
