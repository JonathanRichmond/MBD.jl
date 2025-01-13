"""
CR3BP orbit family wrapper

Author: Jonathan Richmond
C: 1/17/23
U: 1/12/25
"""

import Combinatorics, CSV, DataFrames, LinearAlgebra, StaticArrays
import MBD: CR3BPOrbitFamily

export eigenSort!, exportData, getAlternateIndices, getNumMembers

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
    stabilityIndices::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(orbitFamily.familyMembers))
    alternateStabilityIndices::Vector{Vector{Complex{Float64}}} = Vector{Vector{Complex{Float64}}}(undef, length(orbitFamily.familyMembers))
    for i::Int64 in 1:length(sortedEigenvalues)
        stabilityIndices[i] = 0.5*[LinearAlgebra.norm(sortedEigenvalues[i][1])+1/LinearAlgebra.norm(sortedEigenvalues[i][1]), LinearAlgebra.norm(sortedEigenvalues[i][3])+1/LinearAlgebra.norm(sortedEigenvalues[i][3]), LinearAlgebra.norm(sortedEigenvalues[i][5])+1/LinearAlgebra.norm(sortedEigenvalues[i][5])]
        alternateStabilityIndices[i] = 0.5*[sortedEigenvalues[i][1]+1/sortedEigenvalues[i][1], sortedEigenvalues[i][2]+1/sortedEigenvalues[i][2], sortedEigenvalues[i][3]+1/sortedEigenvalues[i][3]]
    end
    orbitFamily.stabilityIndices = stabilityIndices
    orbitFamily.alternateIndices = alternateStabilityIndices
end

"""
    exportData(orbitFamily, filename)

Export family data to .csv file

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
- `filename::String`: Export file name
"""
function exportData(orbitFamily::CR3BPOrbitFamily, filename::String)
    println("\nExporting family data to file '$filename'...")
    x::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    y::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    z::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    xdot::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    ydot::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    zdot::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    JC::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    p::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    nu::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    tau::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    alpha::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    beta::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    lambda1::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    lambda2::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    lambda3::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    lambda4::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    lambda5::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    lambda6::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    stabilityIndex1::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    stabilityIndex2::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    stabilityIndex3::Vector{Float64} = Vector{Float64}(undef, length(orbitFamily.familyMembers))
    alternateStabilityIndex1::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    alternateStabilityIndex2::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    alternateStabilityIndex3::Vector{Complex{Float64}} = Vector{Complex{Float64}}(undef, length(orbitFamily.familyMembers))
    for m::Int64 in 1:length(orbitFamily.familyMembers)
        x[m] = orbitFamily.familyMembers[m].initialCondition[1]
        y[m] = orbitFamily.familyMembers[m].initialCondition[2]
        z[m] = orbitFamily.familyMembers[m].initialCondition[3]
        xdot[m] = orbitFamily.familyMembers[m].initialCondition[4]
        ydot[m] = orbitFamily.familyMembers[m].initialCondition[5]
        zdot[m] = orbitFamily.familyMembers[m].initialCondition[6]
        JC[m] = orbitFamily.familyMembers[m].JacobiConstant
        p[m] = orbitFamily.familyMembers[m].period
        nu[m] = orbitFamily.familyMembers[m].nu
        tau[m] = orbitFamily.familyMembers[m].tau
        alpha[m] = orbitFamily.familyMembers[m].BrouckeStability[1]
        beta[m] = orbitFamily.familyMembers[m].BrouckeStability[2]
        lambda1[m] = orbitFamily.eigenvalues[m][1]
        lambda2[m] = orbitFamily.eigenvalues[m][2]
        lambda3[m] = orbitFamily.eigenvalues[m][3]
        lambda4[m] = orbitFamily.eigenvalues[m][4]
        lambda5[m] = orbitFamily.eigenvalues[m][5]
        lambda6[m] = orbitFamily.eigenvalues[m][6]
        stabilityIndex1[m] = orbitFamily.stabilityIndices[m][1]
        stabilityIndex2[m] = orbitFamily.stabilityIndices[m][2]
        stabilityIndex3[m] = orbitFamily.stabilityIndices[m][3]
        alternateStabilityIndex1[m] = orbitFamily.alternateIndices[m][1]
        alternateStabilityIndex2[m] = orbitFamily.alternateIndices[m][2]
        alternateStabilityIndex3[m] = orbitFamily.alternateIndices[m][3]
    end
    familyData::DataFrames.DataFrame = DataFrames.DataFrame("x" => x, "y" => y, "z" => z, "xdot" => xdot, "ydot" => ydot, "zdot" => zdot, "JC" => JC, "Period" => p, "Stability Index" => nu, "Time Constant" => tau, "Broucke Stability Parameter 1" => alpha, "Broucke Stability Parameter 2" => beta, "Eigenvalue 1" => lambda1, "Eigenvalue 2" => lambda2, "Eigenvalue 3" => lambda3, "Eigenvalue 4" => lambda4, "Eigenvalue 5" => lambda5, "Eigenvalue 6" => lambda6, "Stability Index 1" => stabilityIndex1, "Stability Index 2" => stabilityIndex2, "Stability Index 3" => stabilityIndex3, "Alternate Stability Index 1" => alternateStabilityIndex1, "Alternate Stability Index 2" => alternateStabilityIndex2, "Alternate Stability Index 3" => alternateStabilityIndex3)
    CSV.write(filename, familyData)
end

"""
    getAlternateIndices(orbitFamily)

Return alternate stability indices for the orbit family

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function getAlternateIndices(orbitFamily::CR3BPOrbitFamily)
    (eigenvalues::Vector{SVector{Complex{Float64}}}, eigenvectors::Vector{SMatrix{Complex{Float64}}}) = eigenSort!(orbitFamily)
    numMembers::Int16 = getNumMembers(orbitFamily)
    standardStabilityIndices::Vector{SVector{Float64}} = [zeros(SVector{3, Float64}) for _ in 1:numMembers]
    alternateStabilityIndices::Vector{SVector{Complex{Float64}}} = [zeros(SVector{3, Complex{Float64}}) for _ in 1:numMembers]
    for i::Int16 in 1:numMembers
        standardStabilityIndices[i] = SVector{3, Float64}(0.5*[LinearAlgebra.norm(orbitFamily.eigenvalues[i][1])+1/LinearAlgebra.norm(orbitFamily.eigenvalues[i][1]), LinearAlgebra.norm(orbitFamily.eigenvalues[i][3])+1/LinearAlgebra.norm(orbitFamily.eigenvalues[i][3]), LinearAlgebra.norm(orbitFamily.eigenvalues[i][5])+1/LinearAlgebra.norm(orbitFamily.eigenvalues[i][5])])
        alternateStabilityIndices[i] = SVector{3, Complex{Float64}}(0.5*[orbitFamily.eigenvalues[i][1]+1/orbitFamily.eigenvalues[i][1], orbitFamily.eigenvalues[i][2]+1/orbitFamily.eigenvalues[i][2], orbitFamily.eigenvalues[i][3]+1/orbitFamily.eigenvalues[i][3]])
    end

    return (standardStabilityIndices, alternateStabilityIndices)
end

"""
    getNumMembers(orbitFamily)

Return number of family members

# Arguments
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function getNumMembers(orbitFamily::CR3BPOrbitFamily)
    return Int16(length(orbitFamily.periods))
end
