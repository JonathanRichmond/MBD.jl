"""
Constraint vector L2 norm convergence check wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 6/9/25
"""

import LinearAlgebra, Logging
import MBD: ConstraintVectorL2NormConvergenceCheck

export isConverged

"""
    isConverged(convergenceCheck, problem)

Return true if problem is converged

# Arguments
- `convergenceCheck::ConstraintVectorL2NormConvergenceCheck`: Constraint vector L2 norm convergence check object
- `problem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function isConverged(convergenceCheck::ConstraintVectorL2NormConvergenceCheck, problem::MBD.BCR4BP12MultipleShooterProblem)
    constraintVec::Vector{Float64} = getConstraintVector!(problem)
    constraintNorm::Float64 = LinearAlgebra.norm(constraintVec)
    maxNorm::Float64 = convergenceCheck.maxVectorNorm

    Logging.@debug "Constraint vector L2 norm: $constraintNorm (threshold: $maxNorm)"

    if constraintNorm <= maxNorm
        Logging.@info "Convergence check passed: norm = $constraintNorm ≤ $maxNorm"

        return true
    else
        Logging.@debug "Convergence check failed: norm = $constraintNorm > $maxNorm"

        return false
    end
end

"""
    isConverged(convergenceCheck, problem)

Return true if problem is converged

# Arguments
- `convergenceCheck::ConstraintVectorL2NormConvergenceCheck`: Constraint vector L2 norm convergence check object
- `problem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function isConverged(convergenceCheck::ConstraintVectorL2NormConvergenceCheck, problem::MBD.CR3BPMultipleShooterProblem)
    constraintVec::Vector{Float64} = getConstraintVector!(problem)
    constraintNorm::Float64 = LinearAlgebra.norm(constraintVec)
    maxNorm::Float64 = convergenceCheck.maxVectorNorm

    Logging.@debug "Constraint vector L2 norm: $constraintNorm (threshold: $maxNorm)"

    if constraintNorm <= maxNorm
        Logging.@info "Convergence check passed: norm = $constraintNorm ≤ $maxNorm"

        return true
    else
        Logging.@debug "Convergence check failed: norm = $constraintNorm > $maxNorm"

        return false
    end
end
