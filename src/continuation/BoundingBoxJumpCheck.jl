"""
Bounding box jump check wrapper

Author: Jonathan Richmond
C: 1/8/23
U: 6/7/25
"""

import Logging
import MBD: BoundingBoxJumpCheck

export addBounds!, checkBounds, isFamilyMember, removeBounds!

"""
    addBounds!(jumpCheck, problem, variable, bounds)

Return bounding box jump check object with updated bounds

# Arguments
- `jumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `problem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `variable::Variable`: Bounded free variable
- `bounds::Matrix{Float64}`: Minimum/maximum values for each free variable
"""
function addBounds!(jumpCheck::BoundingBoxJumpCheck, problem::MBD.CR3BPMultipleShooterProblem, variable::MBD.Variable, bounds::Matrix{Float64})
    index0::Int16 = get(getFreeVariableIndexMap!(problem), variable) do
        err::String = "Variable $(variable.name) is not part of the problem"
        Logging.@error err
        throw(ArgumentError(err))
    end
    numFreeVars::Int64 = getNumFreeVariables(variable)

    Logging.@debug "Adding bounds for variable $(variable.name), indices $index0 to $(index0+numFreeVars-1)"

    checkBounds(jumpCheck, variable, bounds)
    for j in 1:numFreeVars
        minBound::Float64, maxBound::Float64 = bounds[j,1], bounds[j,2]
        if !isnan(minBound) && !isnan(maxBound)
            index::Int16 = index0+j-1
            jumpCheck.variableBounds[index] = copy(bounds[j,:])
            Logging.@debug "Set bounds for index $index: [$minBound, $maxBound]"
        end
    end

    Logging.@debug "Variable $(variable.name) bounds added"
end

"""
    checkBounds(jumpCheck, variable, bounds)

Return error if bounds are invalid

# Arguments
- `jumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `variable::Variable`: Bounded free variable
- `bounds::Matrix{Float64}`: Minimum/maximum values for each free variable
"""
function checkBounds(jumpCheck::BoundingBoxJumpCheck, variable::MBD.Variable, bounds::Matrix{Float64})
    numBounds::Int64 = size(bounds, 1)
    numFreeVars::Int64 = getNumFreeVariables(variable)

    Logging.@debug "Checking bounds for variable $(variable.name): $numBounds bounds vs. $numFreeVars free variables"

    if numBounds != numFreeVars
        err1::String = "Expected $numFreeVars bound rows, found $numBounds"
        Logging.@error err1
        throw(ArgumentError(err1))
    end
    for j in 1:numBounds
        row::Vector{Float64} = bounds[j,:]
        if length(row) != 2
            err2::String = "Row $j has length $(length(row)); expected 2"
            Logging.@error err2
            throw(ArgumentError(err2))
        end
        minBound::Float64, maxBound::Float64 = row[1], row[2]
        if !isnan(minBound) && !isnan(maxBound) && (minBound > maxBound)
            err3::String = "In row $j: minimum bound $minBound > maximum bound $maxBound"
            Logging.@error err3
            throw(ArgumentError(err3))
        end
        Logging.@debug "Row $j bounds OK: [$minBound, $maxBound]"
    end

    Logging.@info "All bounds for variable $(variable.name) passed validation"
end

"""
    isFamilyMember(jumpCheck, data)

Return true if converged solution is family member

# Arguments
- `jumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isFamilyMember(jumpCheck::BoundingBoxJumpCheck, data::MBD.CR3BPContinuationData)
    freeVars::Vector{Float64} = getFreeVariableVector!(data.previousSolution)

    Logging.@debug "Checking if current solution is within family bounding box"

    for (index::Int16, varBounds::Vector{Float64}) in jumpCheck.variableBounds
        value::Float64 = freeVars[index]
        minBound::Float64, maxBound::Float64 = varBounds[1], varBounds[2]
        Logging.@debug "Index $index: value = $value, bounds = [$minBound, $maxBound]"
        if (value < minBound) || (value > maxBound)
            Logging.@info "Bounding box exceeded at index $index: $value ∉ [$(varBounds[1]), $(varBounds[2])]"
            
            return false
        end
    end

    Logging.@info "All free variables are within family bounding box"
    return true
end

"""
    removeBounds!(jumpcheck, problem, variable)

Return bounding box jump check object with updated bounds

# Arguments
- `jumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `problem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `variable::Variable`: Bounded free variable
"""
function removeBounds!(jumpCheck::BoundingBoxJumpCheck, problem::MBD.CR3BPMultipleShooterProblem, variable::MBD.Variable)
    index0::Int16 = get(getFreeVariableIndexMap!(problem), variable) do
        err::String = "Variable $(variable.name) is not part of the problem"
        Logging.@error err
        throw(ArgumentError(err))
    end
    numFreeVars::Int64 = getNumFreeVariables(variable)

    Logging.@debug "Removing bounds for variable $(variable.name), indices $index0 to $(index0+numFreeVars-1)"

    for j in 1:numFreeVars
        index::Int16 = index0+j-1
        if haskey(jumpCheck.variableBounds, index)
            delete!(jumpCheck.variableBounds, index)
            Logging.@debug "Deleted bounds at index $index"
        else
            Logging.@debug "No bounds to delete at index $index"
        end
    end

    Logging.@debug "Variable $(variable.name) bounds removed"
end
