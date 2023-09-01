"""
Minimum norm update generator wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 8/30/23
"""

import LinearAlgebra
import MBD: MinimumNormUpdateGenerator

export canGenerateUpdate, getFullUpdate

"""
    canGenerateUpdate(minimumNormUpdateGenerator, multipleShooterProblem)

Return true if generator can update problem

# Arguments
- `minimumNormUpdateGenerator::MinimumNormUpdateGenerator`: Minimum norm update generator object
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function canGenerateUpdate(minimumNormUpdateGenerator::MinimumNormUpdateGenerator, multipleShooterProblem::MBD.MultipleShooterProblem)
    return getNumberConstraints(multipleShooterProblem) <= getNumberFreeVariables(multipleShooterProblem)
end

"""
    getFullUpdate(minimumNormUpdateGenerator, multipleShooterProblem)

Return update step for free variable vector

# Arguments
- `minimumNormUpdateGenerator::MinimumNormUpdateGenerator`: Minimum norm update generator object
- `multipleShooterProblem::MultipleShooterProblem`: Multiples shooter problem object
"""
function getFullUpdate(minimumNormUpdateGenerator::MinimumNormUpdateGenerator, multipleShooterProblem::MBD.MultipleShooterProblem)
    n_freeVariables::Int64 = getNumberFreeVariables(multipleShooterProblem)
    jacobian::Matrix{Float64} = getJacobian(multipleShooterProblem)
    FX::Vector{Float64} = -1 .*copy(getConstraintVector!(multipleShooterProblem))
    (length(FX) <= n_freeVariables) || throw(ErrorException("Cannot generate update: Number of constraints is greater than number of free variables"))
    solver = LinearAlgebra.qr(jacobian, LinearAlgebra.ColumnNorm())
    X::Vector{Float64} = solver\FX

    return X
end
