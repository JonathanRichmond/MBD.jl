"""
Bifurcation wrapper

Author: Jonathan Richmond
C: 1/18/24
U: 8/15/24
"""

import LinearAlgebra
import MBD: Bifurcation

export getTangentBifurcationStepEigen!, getTangentBifurcationStepSVD!

"""
    getTangentBifurcationStepEigen!(bifurcation)

Return step in tangent bifurcation direction from eigenvectors

# Arguments
- `bifurcation::Bifurcation`: Bifurcation object
"""
function getTangentBifurcationStepEigen!(bifurcation::Bifurcation)
    stepDirection::Vector{Complex{Float64}} = bifurcation.sortedEigenvectors[:,3]
    step::Vector{Float64} = Vector{Float64}(undef, 6)
    for i::Int64 = 1:length(stepDirection)
        (abs(real(stepDirection[i])) >= 1E-4) ? (step[i] = real(stepDirection[i])) : (step[i] = 0.0)
    end
    bifurcation.ICStep = step
end

"""
    getTangentBifurcationStepSVD!(bifurcation)

Return step in tangent bifurcation direction from singular value decomposition

# Arguments
- `bifurcation::Bifurcation`: Bifurcation object
"""
function getTangentBifurcationStepSVD!(bifurcation::Bifurcation)
    DF::Matrix{Float64} = getJacobian(bifurcation.orbit.problem)
    (U::Matrix{Float64}, S::Vector{Float64}, V::Matrix{Float64}) = LinearAlgebra.svd(DF)
    (S[end] < 0.01) || throw(ErrorException("DF matrix has no null space"))
    step::Vector{Float64} = V[:,end]
    (length(step) == getNumberFreeVariables(bifurcation.orbit.problem)) || throw(ErrorException("Step length, $(length(step)), must match number of free variables, $(getNumberFreeVariables(bifurcation.orbit.problem))"))
    for i::Int64 = 1:length(step)
        (abs(step[i]) < 1E-9) && (step[i] = 0.0)
    end
    bifurcation.FVVStep = step
end
