"""
Adaptive step size by element generator wrapper

Author: Jonathan Richmond
C: 1/5/23
U: 6/6/25
"""

import Logging
import MBD: AdaptiveStepSizeByElementGenerator

export updateStepSize!

"""
    updateStepSize!(adaptiveStepSizeByElementGenerator, data)

Return continuation data with updated step size

# Arguments
- `adaptiveStepSizeByElementGenerator::AdaptiveStepSizeByElementGenerator`: Adaptive step size by element generator object
- `data::CR3BPContinuationData`: CR3BP continuation data object
"""
function updateStepSize!(generator::AdaptiveStepSizeByElementGenerator, data::MBD.CR3BPContinuationData)
    step::Float64 = data.currentStepSize
    signFactor::Int64 = (step < 0) ? -1 : 1
    absStep::Float64 = abs(step)

    Logging.@debug "Starting updateStepSize! - Step: $step, Converging: $(data.converging), Iterations: $(data.numIterations)"

    if data.converging
        if data.numIterations < generator.maxIterations
            absStep *= generator.scaleFactor
            absStep = min(absStep, abs(generator.maxStepSize))
            Logging.@debug "Increasing step size to $absStep (scaled up)"
        elseif data.numIterations > generator.minIterations
            absStep /= generator.scaleFactor
            absStep = max(absStep, abs(generator.minStepSize))
            Logging.@debug "Decreasing step size to $absStep (scaled down)"
        else
            Logging.@debug "Step size unchanged: $absStep (iterations within target range)"
        end
        if absStep > abs(generator.maxElementStepSize)
            Logging.@warn "Step size $absStep exceeds maximum element step size $(abs(generator.maxElementStepSize)); clipping"
            absStep = abs(generator.maxElementStepSize)
        end
        data.currentStepSize = signFactor*absStep
    else
        relTol::Float64 = abs((absStep-abs(generator.minStepSize))/generator.minStepSize)
        if relTol < 1E-4
            data.forceEndContinuation = true
            Logging.@info "Terminating continuation: step size $absStep near minimum step size $(abs(generator.minStepSize))"
        else
            absStep = max(abs(generator.minStepSize), absStep/generator.scaleFactor)
            data.currentStepSize = signFactor*absStep
            Logging.@debug "Non-convergent case: reduced step size to $absStep"
        end
    end

    Logging.@debug "Final step size set to: $(data.currentStepSize)"
end
