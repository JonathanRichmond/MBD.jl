"""
CR3BP orbit family wrapper

Author: Jonathan Richmond
C: 1/17/23
U: 9/10/23
"""

import Combinatorics, LinearAlgebra
import MBD: CR3BPOrbitFamily

export eigenSort!

"""
    eigenSort!(orbitFamily)

Return orbit family object with updated eigenvalues/vectors

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function eigenSort!(orbitFamily::CR3BPOrbitFamily)
    eigenvalues::Vector{Vector{Complex{Float64}}} = Vector{Vector{Complex{Float64}}}(undef, length(orbitFamily.familyMembers))
    eigenvectors::Vector{Matrix{Complex{Float64}}} = Vector{Matrix{Complex{Float64}}}(undef, length(orbitFamily.familyMembers))
    for i::Int64 in 1:length(orbitFamily.familyMembers)
        eigenvalues[i] = orbitFamily.familyMembers[i].eigenvalues
        eigenvectors[i] = orbitFamily.familyMembers[i].eigenvectors
    end
    iszero(length(eigenvalues)) && throw(ArgumentError("No eigenvalues to sort"))
    (length(eigenvalues) == length(eigenvectors)) || throw(ArgumentError("There are $(length(eigenvalues)) eigenvalue sets, but there should be $(length(eigenvectors))"))
    println("Sorting $(length(eigenvalues)) sets of eigenvalues/vectors...")
    sortedIndices::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, length(eigenvalues))
    sortedEigenvalues::Vector{Vector{Complex{Float64}}} = Vector{Vector{Complex{Float64}}}(undef, length(eigenvalues))
    sortedEigenvectors::Vector{Matrix{Complex{Float64}}} = Vector{Matrix{Complex{Float64}}}(undef, length(eigenvectors))
    perms::Vector{Vector{Int64}} = collect(Combinatorics.permutations([1, 2, 3, 4, 5, 6]))
    cost::Vector{Float64} = zeros(Float64, length(perms))
    for p::Int64 in 1:length(perms)
        [cost[p] += abs(1.0-eigenvalues[1][perms[p][2*i-1]]*eigenvalues[1][perms[p][2*i]]) for i in 1:3]
    end
    minCostIndex::Int64 = argmin(cost)
    sortedIndices[1] = copy(perms[minCostIndex])
    previousSortedEigenvalues::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, 6)
    previousSortedEigenvectors::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, 6, 6)
    for i::Int64 in 1:6
        previousSortedEigenvalues[i] = eigenvalues[1][sortedIndices[1][i]]
        previousSortedEigenvectors[:,i] = eigenvectors[1][:,sortedIndices[1][i]]
    end
    sortedEigenvalues[1] = copy(previousSortedEigenvalues)
    sortedEigenvectors[1] = copy(previousSortedEigenvectors)
    for i::Int64 in 2:length(eigenvalues)
        dpError::Matrix{Float64} = zeros(Float64, 6, 6)
        for j::Int64 in 1:6
            [dpError[j,k] = abs(1.0-abs(LinearAlgebra.dot(previousSortedEigenvectors[:,j], eigenvectors[i][:,k]))) for k in 1:6]
        end
        cost = zeros(Float64, length(perms))
        for p::Int64 in 1:length(perms)
            for j::Int64 in 1:6
                cost[p] += dpError[j,perms[p][j]]+abs(eigenvalues[i][perms[p][j]]-previousSortedEigenvalues[j])
                any([j == 2, j == 4, j == 6]) && (cost[p] += abs(1.0-eigenvalues[i][perms[p][j-1]]*eigenvalues[i][perms[p][j]]))
            end
        end
        minCostIndex = argmin(cost)
        sortedIndices[i] = copy(perms[minCostIndex])
        for j::Int64 in 1:6
            previousSortedEigenvalues[j] = eigenvalues[i][sortedIndices[i][j]]
            previousSortedEigenvectors[:,j] = eigenvectors[i][:,sortedIndices[i][j]]
        end
        sortedEigenvalues[i] = copy(previousSortedEigenvalues)
        sortedEigenvectors[i] = copy(previousSortedEigenvectors)
    end
    orbitFamily.eigenvalues = sortedEigenvalues
    orbitFamily.eigenvectors = sortedEigenvectors
end
