"""
Variable wrapper

Author: Jonathan Richmond
C: 9/5/22
U: 1/15/25
"""

import MBD: Variable

export getData, getFreeVariableData, getFreeVariableMask, getNumFreeVariables, maskData
export setFreeVariableData!, setFreeVariableMask!

"""
    deepClone(variable)

Return deep copy of variable object

# Arguments
- `variable::Variable`: Variable object
"""
function deepClone(variable::Variable)
    object = Variable(copy(variable.data), copy(variable.freeVariableMask))
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
    return variable.data[variable.freeVariableMask]
end

"""
    getFreeVariableMask(variable)

Return free variable mask

# Arguments
- `variable::Variable`: Variable object
"""
function getFreeVariableMask(variable::Variable)
    return copy(variable.freeVariableMask)
end

"""
    getNumFreeVariables(variable)

Return number of free variables

# Arguments
- `variable::Variable`: Variable object
"""
function getNumFreeVariables(variable::Variable)
    return length(filter(x -> x == true, variable.freeVariableMask))
end

"""
    setFreeVariableData!(variable, freeVariables)

Return variable object with new free variable data

# Arguments
- `variable::Variable`: Variable object
- `freeVariables::Vector{Float64}`: Free variable data
"""
function setFreeVariableData!(variable::Variable, freeVariables::Vector{Float64})
    numFreeVar::Int16 = Int16(getNumFreeVariables(variable))
    (Int16(length(freeVariables)) > numFreeVar) && throw(ArgumentError("Free variable data has $(length(freeVariables)) elements, but should have $numFreeVar"))
    dataIndex::Int16 = Int16(1)
    for maskIndex::Int16 = Int16(1):Int16(length(variable.freeVariableMask))
        if variable.freeVariableMask[maskIndex]
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
    (length(mask) == length(variable.freeVariableMask)) || throw(ArgumentError("Free variable mask has $(length(mask)) elements, but should have $(length(variable.freeVariableMask))"))
    variable.freeVariableMask = mask
end
