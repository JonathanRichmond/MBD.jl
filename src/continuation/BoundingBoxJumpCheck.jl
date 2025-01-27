"""
Bounding box jump check wrapper

Author: Jonathan Richmond
C: 1/8/23
U: 1/26/25
"""

import MBD: BoundingBoxJumpCheck

export addBounds!, checkBounds, isFamilyMember, removeBounds!

"""
    addBounds!(boundingBoxJumpCheck, problem, variable, bounds)

Return bounding box jump check object with updated bounds

# Arguments
- `boundingBoxJumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `problem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `variable::Variable`: Bounded free variable
- `bounds::Matrix{Float64}`: Minimum/maximum values for each free variable
"""
function addBounds!(boundingBoxJumpCheck::BoundingBoxJumpCheck, problem::MBD.CR3BPMultipleShooterProblem, variable::MBD.Variable, bounds::Matrix{Float64})
    index0::Int16 = Int16(get(getFreeVariableIndexMap!(problem), variable, ArgumentError("Variable is not part of problem")))
    checkBounds(boundingBoxJumpCheck, variable, bounds)
    for i::Int16 in Int16(1):Int16(getNumFreeVariables(variable))
        if (!isnan(bounds[i,1]) && !isnan(bounds[i,2]))
            boundingBoxJumpCheck.variableBounds[index0+i-1] = copy(bounds[i,:])
        end
    end
end

"""
    checkBounds(boundingBoxJumpCheck, variable, bounds)

Return error if bounds are invalid

# Arguments
- `boundingBoxJumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `variable::Variable`: Bounded free variable
- `bounds::Matrix{Float64}`: Minimum/maximum values for each free variable
"""
function checkBounds(boundingBoxJumpCheck::BoundingBoxJumpCheck, variable::MBD.Variable, bounds::Matrix{Float64})
    (size(bounds, 1) == getNumFreeVariables(variable)) || throw(ArgumentError("There are $(size(bounds, 1)) boundary entries but there are $(getNumFreeVariables(variable)) free variables"))
    for i::Int16 in Int16(1):Int16(size(bounds, 1))
        (length(bounds[i,:]) == 2) || throw(ArgumentError("There are $(length(bounds[i,:])) boundary values in row $i but there should be 2"))
        if (!isnan(bounds[i,1]) && !isnan(bounds[i,2]) && (bounds[i,1] > bounds[i,2]))
            throw(ArgumentError("Maximum bound must be greater than minimum bound (row $i)"))
        end
    end
end

"""
    isFamilyMember(boundingBoxJumpCheck, data)

Return true if converged solution is family member

# Arguments
- `boundingBoxJumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isFamilyMember(boundingBoxJumpCheck::BoundingBoxJumpCheck, data::MBD.CR3BPContinuationData)
    freeVariables::Vector{Float64} = getFreeVariableVector!(data.previousSolution)
    for (index::Int16, value::Vector{Float64}) in boundingBoxJumpCheck.variableBounds
        freeVariableEntry::Float64 = freeVariables[index]
        if ((freeVariableEntry < value[1]) || (freeVariableEntry > value[2]))
            return false
        end
    end

    return true
end

"""
    removeBounds!(boundingBoxJumpcheck, problem, variable)

Return bounding box jump check object with updated bounds

# Arguments
- `boundingBoxJumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `problem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `variable::Variable`: Bounded free variable
"""
function removeBounds!(boundingBoxJumpCheck::BoundingBoxJumpCheck, problem::MBD.CR3BPMultipleShooterProblem, variable::MBD.Variable)
    index0::Int16 = Int16(get(getFreeVariableIndexMap!(problem), variable, ArgumentError("Variable is not part of problem")))
    [delete!(boundingBoxJumpCheck.variableBounds, index0+i-1) for i in 1:getNumFreeVariables(variable)]
end
