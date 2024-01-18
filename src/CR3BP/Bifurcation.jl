"""
CR3BP bifurcation wrapper

Author: Jonathan Richmond
C: 1/18/24
"""

import MBD: CR3BPBifurcation

export getTangentBifurcationStep

"""
    getTangentBifurcationStep(bifurcation)

Return step in tangent bifurcation direction

# Arguments
- `bifurcation::CR3BPBifurcation`: CR3BP bifurcation object
"""
function getTangentBifurcationStep(bifurcation::CR3BPBifurcation)
    stepDirection::Vector{Complex{Float64}} = bifurcation.sortedEigenvectors[:,3]
    step::Vector{Float64} = Vector{Float64}(undef, 6)
    for i::Int64 = 1:length(stepDirection)
        (real(stepDirection[i]) >= 1E-4) ? (step[i] = real(stepDirection[i])) : (step[i] = 0.0)
    end

    return step
end
