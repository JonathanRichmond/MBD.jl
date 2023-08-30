"""
Least squares update generator wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 8/30/23
"""

import LinearAlgebra
import MBD: LeastSquaresUpdateGenerator

export canGenerateUpdate, getFullUpdate

"""
    canGenerateUpdate(leastSquaresUpdateGenerator, multipleShooterProblem)

Return true if generator can update problem

# Arguments
- `leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator`: Least squares update generator object
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function canGenerateUpdate(leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator, multipleShooterProblem::MBD.MultipleShooterProblem)
    return getNumberConstraints(multipleShooterProblem) > getNumberFreeVariables(multipleShooterProblem)
end

"""
    getFullUpdate(leastSquaresUpdateGenerator, multipleShooterProblem)

Return update step for free variable vector

# Arguments
- `leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator`: Least squares update generator object
- `multipleShooterProblem::MultipleShooterProblem`: Multiples shooter problem object
"""
function getFullUpdate(leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator, multipleShooterProblem::MBD.MultipleShooterProblem)
    n_freeVariables::Int64 = getNumberFreeVariables(multipleShooterProblem)
    jacobian::Matrix{Float64} = getJacobian(multipleShooterProblem)
    FX::Vector{Float64} = -1 .*copy(getConstraintVector!(multipleShooterProblem))
    (length(FX) > n_freeVariables) || throw(ErrorException("Cannot generate update: Number of constraints is not greater than number of free variables"))
    solver = LinearAlgebra.qr(jacobian'*jacobian)
    X::Vector{Float64} = solver\(jacobian'*FX)

    return X
end