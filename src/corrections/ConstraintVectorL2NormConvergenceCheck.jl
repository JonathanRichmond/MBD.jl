"""
Constraint vector L2 norm convergence check wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 8/30/23
"""

import LinearAlgebra
import MBD: ConstraintVectorL2NormConvergenceCheck

export isConverged

"""
    isConverged(constraintVectorL2NormConvergenceCheck, mulitpleShooterProblem)

Return true if problem is converged

# Arguments
- `constraintVectorL2NormConvergenceCheck::ConstraintVectorL2NormConvergenceCheck`: Constraint vector L2 norm convergence check object
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function isConverged(constraintVectorL2NormConvergenceCheck::ConstraintVectorL2NormConvergenceCheck, multipleShooterProblem::MBD.MultipleShooterProblem)
    return LinearAlgebra.norm(getConstraintVector!(multipleShooterProblem)) <= constraintVectorL2NormConvergenceCheck.maxVectorNorm
end