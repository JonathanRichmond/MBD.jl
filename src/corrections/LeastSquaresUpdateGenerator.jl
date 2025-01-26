"""
Least squares update generator wrapper

Author: Jonathan Richmond
C: 9/9/22
U: 1/16/25
"""

import LinearAlgebra, StaticArrays
import MBD: LeastSquaresUpdateGenerator

export canGenerateUpdate, getFullUpdate

"""
    canGenerateUpdate(leastSquaresUpdateGenerator, multipleShooterProblem)

Return true if generator can update problem

# Arguments
- `leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator`: Least squares update generator object
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function canGenerateUpdate(leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator, multipleShooterProblem::MBD.CR3BPMultipleShooterProblem)
    return getNumConstraints(multipleShooterProblem) > getNumFreeVariables!(multipleShooterProblem)
end

"""
    getFullUpdate(leastSquaresUpdateGenerator, multipleShooterProblem)

Return update step for free variable vector

# Arguments
- `leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator`: Least squares update generator object
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiples shooter problem object
"""
function getFullUpdate(leastSquaresUpdateGenerator::LeastSquaresUpdateGenerator, multipleShooterProblem::MBD.CR3BPMultipleShooterProblem)
    numConstraints::Int64 = getNumConstraints(multipleShooterProblem)
    numFreeVariables::Int64 = getNumFreeVariables!(multipleShooterProblem)
    jacobian::Matrix{Float64} = getJacobian!(multipleShooterProblem)
    FX::StaticArrays.SVector{numConstraints, Float64} = -1 .*getConstraintVector!(multipleShooterProblem)
    (length(FX) > numFreeVariables) || throw(ErrorException("Cannot generate update: Number of constraints is not greater than number of free variables"))
    solver = LinearAlgebra.qr(jacobian'*jacobian)
    X::Vector{Float64} = solver\(jacobian'*FX)

    return X
end
