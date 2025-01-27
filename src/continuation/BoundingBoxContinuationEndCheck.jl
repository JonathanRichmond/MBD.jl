"""
Bounding box continuation end check wrapper

Author: Jonathan Richmond
C: 1/9/23
U: 1/26/25
"""

import MBD: BoundingBoxContinuationEndCheck

export checkBounds, isContinuationDone

"""
    checkBounds(boundingBoxContinuationEndCheck, variable)

Return error if bounds are invalid

# Arguments
- `boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `variable::Variable`: Bounded free variable
"""
function checkBounds(boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck, variable::MBD.Variable)
    (size(boundingBoxContinuationEndCheck.paramBounds, 1) == getNumFreeVariables(variable)) || throw(ArgumentError("There are $(size(boundingBoxContinuationEndCheck.paramBounds, 1)) boundary entries but there are $(getNumFreeVariables(variable)) free variables"))
    for i::Int16 in Int16(1):Int16(size(boundingBoxContinuationEndCheck.paramBounds, 1))
        (length(boundingBoxContinuationEndCheck.paramBounds[i,:]) == 2) || throw(ArgumentError("There are $(length(boundingBoxContinuationEndCheck.paramBounds[i,:])) boundary values in row $i but there should be 2"))
        if (!isnan(boundingBoxContinuationEndCheck.paramBounds[i,1]) && !isnan(boundingBoxContinuationEndCheck.paramBounds[i,2]) && (boundingBoxContinuationEndCheck.paramBounds[i,1] > boundingBoxContinuationEndCheck.paramBounds[i,2]))
            throw(ArgumentError("Maximum bound must be greater than minimum bound (row $i)"))
        end
    end
end

"""
    isContinuationDone(boundingBoxContinuationEndCheck, data)

Return true if continuation is done

# Arguments
- `boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck`: Bounding box continuation end check object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function isContinuationDone(boundingBoxContinuationEndCheck::BoundingBoxContinuationEndCheck, data::MBD.CR3BPContinuationData)
    for (index1::MBD.Variable, value1::Int16) in getFreeVariableIndexMap!(data.previousSolution)
        if index1.name == boundingBoxContinuationEndCheck.paramName
            checkBounds(boundingBoxContinuationEndCheck, index1)
            for i::Int16 in Int16(1):Int16(getNumFreeVariables(index1))
                if (!isnan(boundingBoxContinuationEndCheck.paramBounds[i,1]) && !isnan(boundingBoxContinuationEndCheck.paramBounds[i,2]))
                    boundingBoxContinuationEndCheck.variableBounds[value1+i-1] = copy(boundingBoxContinuationEndCheck.paramBounds[i,:])
                end
            end
            freeVariableVector::Vector{Float64} = getFreeVariableVector!(data.previousSolution)
            for (index2::Int16, value2::Vector{Float64}) in boundingBoxContinuationEndCheck.variableBounds
                freeVariable::Float64 = freeVariableVector[index2]
                if ((freeVariable < value2[1]) || (freeVariable > value2[2]))
                    println("Continuation bounding box reached!")
                    return true
                end
            end
            [delete!(boundingBoxContinuationEndCheck.variableBounds, value1+i-1) for i in 1:getNumFreeVariables(index1)]
        end
    end

    return false
end
