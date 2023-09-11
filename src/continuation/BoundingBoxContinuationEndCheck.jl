"""
Bounding box continuation end check wrapper

Author: Jonathan Richmond
C: 1/9/23
U: 9/10/23
"""

import MBD: BoundingBoxContinuationEndCheck

export checkBounds, isContinuationDone

"""
    checkBounds(boundingBoxContinuationEndCheck, variable, bounds)

Return error if bounds are invalid

# Arguments
- `boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `variable::Variable`: Bounded free variable
- `bounds::Matrix{Float64}`: Minimum/maximum values for each free variable
"""
function checkBounds(boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck, variable::MBD.Variable, bounds::Matrix{Float64})
    (size(bounds, 1) == getNumberFreeVariables(variable)) || throw(ArgumentError("There are $(size(bounds, 1)) boundary entries but there are $(getNumberFreeVariables(variable)) free variables"))
    for i::Int64 in 1:size(bounds, 1)
        (length(bounds[i,:]) == 2) || throw(ArgumentError("There are $(length(bounds[i,:])) boundary values in row $i but there should be 2"))
        if (!isnan(bounds[i,1]) && !isnan(bounds[i,2]) && (bounds[i,1] > bounds[i,2]))
            throw(ArgumentError("Maximum bound must be greater than minimum bound (row $i)"))
        end
    end
end

"""
    isContinuationDone(boundingBoxContinuationEndCheck, data)

Return true if continuation is done

# Arguments
- `boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `data::ContinuationData`: Continuation data object
"""
function isContinuationDone(boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck, data::MBD.ContinuationData)
    for (index1::MBD.Variable, value1::Int64) in getFreeVariableIndexMap!(data.previousSolution)
        if index1.name == boundingBoxContinuationEndCheck.paramName
            checkBounds(boundingBoxContinuationEndCheck, index1, boundingBoxContinuationEndCheck.paramBounds)
            for i::Int64 in 1:getNumberFreeVariables(index1)
                if (!isnan(boundingBoxContinuationEndCheck.paramBounds[i,1]) && !isnan(boundingBoxContinuationEndCheck.paramBounds[i,2]))
                    boundingBoxContinuationEndCheck.variableBounds[value1+i-1] = copy(boundingBoxContinuationEndCheck.paramBounds[i,:])
                end
            end
            freeVariableVector::Vector{Float64} = getFreeVariableVector!(data.previousSolution)
            for (index2::Int64, value2::Vector{Float64}) in boundingBoxContinuationEndCheck.variableBounds
                freeVariable::Float64 = freeVariableVector[index2]
                if ((freeVariable < value2[1]) || (freeVariable > value2[2]))
                    println("Continuation bounding box reached!")
                    return true
                end
            end
            [delete!(boundingBoxContinuationEndCheck.variableBounds, value1+i-1) for i in 1:getNumberFreeVariables(index1)]
        end
    end

    return false
end
