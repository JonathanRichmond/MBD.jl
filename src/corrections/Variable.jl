"""
Variable wrapper

Author: Jonathan Richmond
C: 9/5/22
U: 8/5/23
"""

import MBD: Variable

export getData, getFreeVariableData, getFreeVariableMask
export getNumberFreeVariables, maskData

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
