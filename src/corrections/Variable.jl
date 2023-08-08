"""
Variable wrapper

Author: Jonathan Richmond
C: 9/5/22
U: 8/5/23
"""

import MBD: Variable

export getData, getFreeVariableData, getFreeVariableMask
export getNumberFreeVariables, maskData, setFreeVariableData!
export setFreeVariableMask!

"""
    deepClone(variable)

Return deep copy of variable object

# Arguments
- `variable::Variable`: Variable object
"""
function deepClone(variable::Variable)
    object = Variable(copy(variable.data), copy(variable.freeVarMask))
    object.name = variable.name

    return object
end

"""
    getData(variable)

Return data

# Arguments
- `variable::Variable`: Variable object
"""
function getData(variable::Variable)
    return copy(variable.data)
end

"""
    getFreeVariableData(variable)

Return free variable data

# Arguments
- `variable::Variable`: Variable object
"""
function getFreeVariableData(variable::Variable)
    return variable.data[variable.freeVarMask]
end

"""
    getFreeVariableMask(variable)

Return free variable mask

# Arguments
- `variable::Variable`: Variable object
"""
function getFreeVariableMask(variable::Variable)
    return copy(variable.freeVarMask)
end

"""
    getNumberFreeVariables(variable)

Return number of free variables

# Arguments
- `variable::Variable`: Variable object
"""
function getNumberFreeVariables(variable::Variable)
    return length(filter(x -> x == true, variable.freeVarMask))
end

"""
    setFreeVariableData!(variable, freeVariables)

Return variable object with new free variable data

# Arguments
- `variable::Variable`: Variable object
- `freeVariables::Vector{Float64}`: Free variable data
"""
function setFreeVariableData!(variable::Variable, freeVariables::Vector{Float64})
    n_freeVar::Int64 = getNumberFreeVariables(variable)
    (length(freeVariables) > n_freeVar) && throw(ArgumentError("Free variable data has $(length(freeVariables)) elements, but should have $n_freeVar"))
    dataIndex::Int64 = 1
    for maskIndex::Int64 = 1:length(variable.freeVarMask)
        if variable.freeVarMask[maskIndex]
            variable.data[maskIndex] = freeVariables[dataIndex]
            dataIndex += 1
        end
    end
end

"""
    setFreeVariableMask!(variable, mask)

Return variable object with new free variable mask

# Arguments
- `variable::Variable`: Variable object
- `mask::Vector{Bool}`: Free variable mask
"""
function setFreeVariableMask!(variable::Variable, mask::Vector{Bool})
    (length(mask) == length(variable.freeVarMask)) || throw(ArgumentError("Free variable mask has $(length(mask)) elements, but should have $(length(variable.freeVarMask))"))
    variable.freeVarMask = copy(mask)
end
