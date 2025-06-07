"""
Bounding box continuation end check wrapper

Author: Jonathan Richmond
C: 1/9/23
U: 6/7/25
"""

import Logging
import MBD: BoundingBoxContinuationEndCheck

export checkBounds, isContinuationDone

"""
    checkBounds(boundsCheck, variable)

Return error if bounds are invalid

# Arguments
- `boundsCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `variable::Variable`: Bounded free variable
"""
function checkBounds(boundsCheck::BoundingBoxContinuationEndCheck, variable::MBD.Variable)
    bounds::Matrix{Float64} = boundsCheck.paramBounds
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
    isContinuationDone(boundsCheck, data)

Return true if continuation is done

# Arguments
- `boundsCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isContinuationDone(boundsCheck::BoundingBoxContinuationEndCheck, data::MBD.CR3BPContinuationData)
    Logging.@debug "Checking if continuation bounding box is reached"

    for (var::MBD.Variable, index0::Int16) in getFreeVariableIndexMap!(data.previousSolution)
        if var.name == boundsCheck.paramName
            Logging.@debug "Checking bounds for variable $(var.name)"
            checkBounds(boundsCheck, var)
            numFreeVars::Int64 = getNumFreeVariables(var)
            bounds::Matrix{Float64} = boundsCheck.paramBounds
            for j in 1:numFreeVars
                minBound::Float64, maxBound::Float64 = bounds[j,1], bounds[j,2]
                if !isnan(minBound) && !isnan(maxBound)
                    boundsCheck.variableBounds[Int16(index0+j-1)] = copy(bounds[j,:])
                    Logging.@debug "Set bounds for index $(index0+j-1): [$minBound, $maxBound]"
                end
            end
            freeVars::Vector{Float64} = getFreeVariableVector!(data.previousSolution)
            for (index::Int16, varBounds::Vector{Float64}) in boundsCheck.variableBounds
                value::Float64 = freeVars[index]
                minBound::Float64, maxBound::Float64 = varBounds[1], varBounds[2]
                Logging.@debug "Index $index: value = $value, bounds = [$minBound, $maxBound]"
                if (value < minBound) || (value > maxBound)
                    Logging.@info "Continuation bounding box reached at index $index: $value ∉ [$(varBounds[1]), $(varBounds[2])]"
                    println("Continuation bounding box reached!")
                    
                    return true
                end
            end
            [delete!(boundsCheck.variableBounds, index0+j-1) for j in 1:numFreeVars]
        end
    end

    Logging.@debug "Continuation bounding box not yet reached"
    return false
end
