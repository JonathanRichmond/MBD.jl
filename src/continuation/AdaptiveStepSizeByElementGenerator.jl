"""
Adaptive step size by element generator wrapper

Author: Jonathan Richmond
C: 1/5/23
U: 9/6/23
"""

import MBD: AdaptiveStepSizeByElementGenerator

export updateStepSize!

"""
    updateStepSize!(adaptiveStepSizeByElementGenerator, data)

Return adaptive step size by element generator with updated step size

# Arguments
- `adaptiveStepSizeByElementGenerator::AdaptiveStepSizeByElementGenerator`: Adaptive step size by element generator object
- `data::ContinuationData`: Continuation data object
"""
function updateStepSize!(adaptiveStepSizeByElementGenerator::AdaptiveStepSizeByElementGenerator, data::MBD.ContinuationData)
    if data.converging
        signFactor = (data.currentStepSize < 0) ? -1 : 1
        if data.numIterations < adaptiveStepSizeByElementGenerator.maxIterations
            data.currentStepSize = signFactor*min(abs(adaptiveStepSizeByElementGenerator.maxStepSize), abs(data.currentStepSize*adaptiveStepSizeByElementGenerator.scaleFactor))
        elseif data.numIterations > adaptiveStepSizeByElementGenerator.minIterations
            data.currentStepSize = signFactor*max(abs(adaptiveStepSizeByElementGenerator.minStepSize), abs(data.currentStepSize/adaptiveStepSizeByElementGenerator.scaleFactor))
        end
        (abs(data.currentStepSize) > abs(adaptiveStepSizeByElementGenerator.maxElementStepSize)) && (data.currentStepSize = adaptiveStepSizeByElementGenerator.maxElementStepSize)
    else
        if abs(data.currentStepSize-adaptiveStepSizeByElementGenerator.minStepSize)/abs(adaptiveStepSizeByElementGenerator.minStepSize) < 1E-4
            data.forceEndContinuation = true
        else
            data.currentStepSize = signFactor*max(abs(adaptiveStepSizeByElementGenerator.minStepSize), abs(data.currentStepSize/adaptiveStepSizeByElementGenerator.scaleFactor))
        end
    end
end
