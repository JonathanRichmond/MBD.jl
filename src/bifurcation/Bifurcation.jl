"""
Bifurcation wrapper

Author: Jonathan Richmond
C: 1/18/24
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
    stepDirection::Vector{Float64} = V[:,end]
    step::Vector{Float64} = Vector{Float64}(undef, getNumberFreeVariables(bifurcation.orbit.problem))
    (length(stepDirection) == length(step)) || throw(ErrorException("Step length, $(length(stepDirection)), must match number of free variables, $(length(step))"))
    for i::Int64 = 1:length(stepDirection)
        (stepDirection[i] < 1E-6) && (stepDirection[i] = 0.0)
    end
    bifurcation.FVVStep = step
end
