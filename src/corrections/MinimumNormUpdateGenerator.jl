"""
Minimum norm update generator wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 1/16/25
"""

import LinearAlgebra, StaticArrays
import MBD: MinimumNormUpdateGenerator

export canGenerateUpdate, getFullUpdate

"""
    canGenerateUpdate(minimumNormUpdateGenerator, multipleShooterProblem)

Return true if generator can update problem

# Arguments
- `minimumNormUpdateGenerator::MinimumNormUpdateGenerator`: Minimum norm update generator object
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function canGenerateUpdate(minimumNormUpdateGenerator::MinimumNormUpdateGenerator, multipleShooterProblem::MBD.CR3BPMultipleShooterProblem)
    return getNumConstraints(multipleShooterProblem) <= getNumFreeVariables!(multipleShooterProblem)
end

"""
    getFullUpdate(minimumNormUpdateGenerator, multipleShooterProblem)

Return update step for free variable vector

# Arguments
- `minimumNormUpdateGenerator::MinimumNormUpdateGenerator`: Minimum norm update generator object
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getFullUpdate(minimumNormUpdateGenerator::MinimumNormUpdateGenerator, multipleShooterProblem::MBD.CR3BPMultipleShooterProblem)
    numConstraints::Int64 = getNumConstraints(multipleShooterProblem)
    numFreeVariables::Int64 = getNumFreeVariables!(multipleShooterProblem)
    jacobian::Matrix{Float64} = getJacobian!(multipleShooterProblem)
    FX::StaticArrays.SVector{numConstraints, Float64} = -1 .*getConstraintVector!(multipleShooterProblem)
    (length(FX) <= numFreeVariables) || throw(ErrorException("Cannot generate update: Number of constraints is greater than number of free variables"))
    solver = LinearAlgebra.qr(jacobian, LinearAlgebra.ColumnNorm())
    X::Vector{Float64} = solver\FX

    return X
end
