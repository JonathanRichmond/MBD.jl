"""
CR3BP orbit family wrapper

Author: Jonathan Richmond
C: 1/17/23
U: 6/3/25
"""

import Combinatorics, LinearAlgebra, StaticArrays
import MBD: CR3BPOrbitFamily

export eigenSort!, getAlternateIndices, getMember, getNumMembers

"""
    eigenSort!(orbitFamily)

Return orbit family object with updated eigenvalues/vectors

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function eigenSort!(orbitFamily::CR3BPOrbitFamily)
    nMembers::Int16 = Int16(getNumMembers(orbitFamily))
    iszero(nMembers) && throw(ArgumentError("No eigenvalues to sort"))
    eigenvalues::Vector{Vector{Complex{Float64}}} = Vector{Vector{Complex{Float64}}}(undef, nMembers)
    eigenvectors::Vector{Matrix{Complex{Float64}}} = Vector{Matrix{Complex{Float64}}}(undef, nMembers)
    for o::Int16 in Int16(1):nMembers
        orbit = MBD.CR3BPPeriodicOrbit(orbitFamily.dynamicsModel, orbitFamily.initialConditions[o], orbitFamily.periods[o], orbitFamily.monodromies[o])
        (eigenvalues[o], eigenvectors[o]) = getEigenData(orbit)
    end
    println("Sorting $nMembers sets of eigenvalues/vectors...")
    sortedIndices::Vector{StaticArrays.MVector{6, Int64}} = Vector{StaticArrays.MVector{6, Int64}}(undef, nMembers)
    sortedEigenvalues::Vector{Vector{Complex{Float64}}} = Vector{Vector{Complex{Float64}}}(undef, nMembers)
    sortedEigenvectors::Vector{Matrix{Complex{Float64}}} = Vector{Matrix{Complex{Float64}}}(undef, nMembers)
    perms::Vector{Vector{Int64}} = collect(Combinatorics.permutations([1, 2, 3, 4, 5, 6]))
    nPerms::Int16 = Int16(length(perms))
    cost::Vector{Float64} = zeros(Float64, nPerms)
    for p::Int16 in Int16(1):nPerms
        [cost[p] += abs(1.0-eigenvalues[1][perms[p][2*i-1]]*eigenvalues[1][perms[p][2*i]]) for i in 1:3]
    end
    minCostIndex::Int16 = Int16(argmin(cost))
    sortedIndices[1] = StaticArrays.MVector{6, Int64}(perms[minCostIndex])
    previousSortedEigenvalues::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, 6)
    previousSortedEigenvectors::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, 6, 6)
    for i::Int16 in Int16(1):Int16(6)
        previousSortedEigenvalues[i] = eigenvalues[1][sortedIndices[1][i]]
        previousSortedEigenvectors[:,i] = eigenvectors[1][:,sortedIndices[1][i]]
    end
    sortedEigenvalues[1] = copy(previousSortedEigenvalues)
    sortedEigenvectors[1] = copy(previousSortedEigenvectors)
    for i::Int16 in Int16(2):nMembers
        dpError::Matrix{Float64} = zeros(Float64, (6,6))
        for j::Int16 in Int16(1):Int16(6)
            [dpError[j,k] = abs(1.0-abs(LinearAlgebra.dot(previousSortedEigenvectors[:,j], eigenvectors[i][:,k]))) for k in 1:6]
        end
        cost = zeros(Float64, nPerms)
        for p::Int16 in Int16(1):nPerms
            for j::Int16 in Int16(1):Int16(6)
                cost[p] += dpError[j,perms[p][j]]+abs(eigenvalues[i][perms[p][j]]-previousSortedEigenvalues[j])
                any([j == 2, j == 4, j == 6]) && (cost[p] += abs(1.0-eigenvalues[i][perms[p][j-1]]*eigenvalues[i][perms[p][j]]))
            end
        end
        minCostIndex = Int16(argmin(cost))
        sortedIndices[i] = StaticArrays.MVector{6, Int64}(perms[minCostIndex])
        for j::Int16 in Int16(1):Int16(6)
            previousSortedEigenvalues[j] = eigenvalues[i][sortedIndices[i][j]]
            previousSortedEigenvectors[:,j] = eigenvectors[i][:,sortedIndices[i][j]]
        end
        sortedEigenvalues[i] = copy(previousSortedEigenvalues)
        sortedEigenvectors[i] = copy(previousSortedEigenvectors)
    end
    orbitFamily.eigenvalues = sortedEigenvalues
    orbitFamily.eigenvectors = sortedEigenvectors
    orbitFamily.hasBeenSorted = true
end

"""
    getAlternateIndices(orbitFamily)

Return alternate stability indices for the orbit family

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function getAlternateIndices(orbitFamily::CR3BPOrbitFamily)
    orbitFamily.hasBeenSorted || eigenSort!(orbitFamily)
    numMembers::Int16 = getNumMembers(orbitFamily)
    standardStabilityIndices::Vector{Vector{Float64}} = [zeros(Float64, 3) for _ in 1:numMembers]
    alternateStabilityIndices::Vector{Vector{Complex{Float64}}} = [zeros(Complex{Float64}, 3) for _ in 1:numMembers]
    for i::Int16 in Int16(1):numMembers
        standardStabilityIndices[i] = 0.5*[LinearAlgebra.norm(orbitFamily.eigenvalues[i][1])+1/LinearAlgebra.norm(orbitFamily.eigenvalues[i][1]), LinearAlgebra.norm(orbitFamily.eigenvalues[i][3])+1/LinearAlgebra.norm(orbitFamily.eigenvalues[i][3]), LinearAlgebra.norm(orbitFamily.eigenvalues[i][5])+1/LinearAlgebra.norm(orbitFamily.eigenvalues[i][5])]
        alternateStabilityIndices[i] = 0.5*[orbitFamily.eigenvalues[i][1]+1/orbitFamily.eigenvalues[i][1], orbitFamily.eigenvalues[i][3]+1/orbitFamily.eigenvalues[i][3], orbitFamily.eigenvalues[i][5]+1/orbitFamily.eigenvalues[i][5]]
    end

    return (standardStabilityIndices, alternateStabilityIndices)
end

"""
    getMember(orbitFamily, orbit)

Return periodic orbit member object

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
- `orbit::Int64`: Orbit identifier
"""
function getMember(orbitFamily::CR3BPOrbitFamily, orbit::Int64)
    return MBD.CR3BPPeriodicOrbit(orbitFamily.dynamicsModel, orbitFamily.initialConditions[orbit], orbitFamily.periods[orbit], orbitFamily.monodromies[orbit])
end

"""
    getNumMembers(orbitFamily)

Return number of family members

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function getNumMembers(orbitFamily::CR3BPOrbitFamily)
    return length(orbitFamily.periods)
end
