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

    Logging.@debug "Checking parameter bounds: $numBounds bounds vs. $numFreeVars free variables"

    if numBounds != numFreeVars
        err::String = "Expected $numFreeVars bound rows, found $numBounds"
        Logging.@error err
        throw(ArgumentError(err))
    end
    for j in 1:numBounds
        row::Vector{Float64} = bounds[j,:]
        if length(row) != 2
            err::String = "Row $j has length $(length(row)); expected 2"
            Logging.@error err
            throw(ArgumentError(err))
        end
        minBound::Float64, maxBound::Float64 = row[1], row[2]
        if !isnan(min_bound) && !isnan(maxBound) && (minBound > maxBound)
            err::String = "In row $j: minimum bound $minBound > maximum bound $maxBound"
            Logging.@error err
            throw(ArgumentError(err))
        end
        Logging.@debug "Row $j bounds OK: [$minBound, $maxBound]"
    end

    Logging.@info "All parameter bounds passed validation"
end

"""
    isContinuationDone(boundsCheck, data)

Return true if continuation is done

# Arguments
- `boundsCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isContinuationDone(boundsCheck::BoundingBoxContinuationEndCheck, data::MBD.CR3BPContinuationData)
    Logging.@debug "Bounding box continuation end check"

    for (var::MBD.Variable, index::Int16) in getFreeVariableIndexMap!(data.previousSolution)
        if var.name == boundsCheck.paramName
            Logging.@debug "Checking bounds for variable $(var.name)"
            checkBounds(boundsCheck, var)
            numFreeVars::Int64 = getNumFreeVariables(var)
            bounds::Vector{Float64} = boundsCheck.paramBounds
            for j in 1:numFreeVars
                minBounds::Float64, maxBounds::Float64 = bounds[j,1], bounds[j,2]
                if !isnan(minBounds) && !isnan(maxBounds)
                    boundsCheck.variableBounds[Int16(index+j-1)] = copy(bounds[j,:])
                    Logging.@debug "Set bounds for index $(index+j-1): [$minBounds, $maxBounds]"
                end
            end
            freeVars::Vector{Float64} = getFreeVariableVector!(data.previousSolution)
            for (index2::Int16, bounds2::Vector{Float64}) in boundsCheck.variableBounds
                value::Float64 = freeVars[index2]
                if (value < bounds2[1]) || (value > bounds2[2])
                    Logging.@info "Continuation bounding box reached at index $index2: $value ∉ [$(bounds2[1]), $(bounds2[2])]"
                    println("Continuation bounding box reached!")
                    
                    return true
                end
            end
            [delete!(boundsCheck.variableBounds, index+j-1) for j in 1:numFreeVars]
        end
    end

    return false
end
