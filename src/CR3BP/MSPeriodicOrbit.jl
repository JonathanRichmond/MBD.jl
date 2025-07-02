"""
CR3BP multiple shooter periodic orbit wrapper

Author: Jonathan Richmond
C: 7/2/25
"""

import DifferentialEquations, LinearAlgebra, SparseArrays, StaticArrays
import MBD: CR3BPMSPeriodicOrbit

export getApproxEigenData, getBrouckeStability, getEigenData, getJacobiConstant, getStabilityIndex
export getTimeConstant

"""
    getApproxEigenData(periodicOrbit; clusterTol, complexTol)

Return approximate orbit eigenvalues and -vectors

# Arguments
- `periodicOrbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
- `clusterTol::Float64`: Clustering relative tolerance (default = 1E-5)
- `complexTol::Float64`: Complex number relative tolerance (default = 1E-5)
"""
function getApproxEigenData(periodicOrbit::CR3BPMSPeriodicOrbit, clusterTol::Float64 = 1E-5, complexTol::Float64 = 1E-5)
    propagator = MBD.Propagator(equationType = MBD.STM)
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    nSegs::Int64 = length(periodicOrbit.nodeEpochs)-1
    renormalizeEvent::DifferentialEquations.DiscreteCallback = DifferentialEquations.PeriodicCallback(MBD.renormalize!, 0.1)
    STMs::Vector{Matrix{Float64}} = []
    for s::Int64 = 1:nSegs
        Rs::Vector{Matrix{Float64}} = []
        propSegment::MBD.CR3BPArc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(periodicOrbit.dynamicsModel, periodicOrbit.nodeStates[s], MBD.STM), [periodicOrbit.nodeEpochs[s] periodicOrbit.nodeEpochs[s+1]], targeter.dynamicsModel, [targeter.dynamicsModel, Rs])
        endState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(propSegment, -1))
        STM::Matrix{Float64} = reshape(endState[7:42], (6,6))
        for R::Matrix{Float64} in reverse(Rs)
            STM *= R
        end
        push!(STMs, STM)
    end
    rows::Vector{Int64} = []
    cols::Vector{Int64} = []
    vals::Vector{Float64} = []
    r_offset::Int64 = 0
    c_offset::Int64 = 0
    for S::StaticArrays.SMatrix{6, 6, Float64} in STMs
        for i::Int64 = 1:6, j::Int64 = 1:6
            push!(rows, r_offset+i)
            push!(cols, c_offset+j)
            push!(vals, S[i,j])
        end
        r_offset += 6
        c_offset += 6
    end
    Phi::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.sparse(rows, cols, vals)
    offI::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(2*nSegs*6, 2*nSegs*6)
    for i::Int64 = 1:2*nSegs-1
        offI[(6*(i-1)+1):(6*i),(6*i+1):(6*(i+1))] = SparseArrays.sparse(LinearAlgebra.I, 6, 6)
    end
    offI[(2*nSegs*6-5):(2*nSegs*6),1:6] = SparseArrays.sparse(LinearAlgebra.I, 6, 6)
    E::LinearAlgebra.Eigen = LinearAlgebra.eigen(offI\Matrix(Phi))
    Lambda::Vector{Complex{Float64}} = E.values.^(2*nSegs)
    indices::Vector{Int64} = sortperm(Lambda, by = x -> abs(real(x)))
    sortedLambda::Vector{Complex{Float64}} = Lambda[indices]
    sortedV::Matrix{Complex{Float64}} = E.vectors[1:6,indices]
    clusters::Vector{Vector{Tuple{Complex{Float64}, Vector{Complex{Float64}}}}} = [[(sortedLambda[1], sortedV[:,1])]]
    for e::Int64 in 2:(2*nSegs*6)
        lambda::Complex{Float64} = sortedLambda[e]
        v::Vector{Complex{Float64}} = sortedV[:,e]
        added::Bool = false
        for c::Vector{Tuple{Complex{Float64}, Vector{Complex{Float64}}}} in clusters
            lambdas::Vector{Complex{Float64}} = [x[1] for x in c]
            isConjugate::Bool = (any(x -> isapprox(lambda, conj(x); atol = 1E-8), lambdas) && (abs(imag(lambda)) > complexTol*max(abs(real(lambda)), 1E-12)))
            if (abs(real(lambda)-real(lambdas[end]))/max(abs(real(lambda)), abs(real(lambdas[end])), 1E-12) < clusterTol) && !isConjugate
                push!(c, (lambda, v))
                added = true
                break
            end
        end
        !added && push!(clusters, [(lambda, v)])
    end
    Lambda_avg::Vector{Complex{Float64}} = Vector{Float64}(undef, length(clusters))
    V_avg::Matrix{Complex{Float64}} = Matrix{Float64}(undef, 6, length(clusters))
    for c::Int64 = 1:length(clusters)
        lambdas::Vector{Complex{Float64}} = [x[1] for x in clusters[c]]
        vs::Vector{Vector{Complex{Float64}}} = [x[2] for x in clusters[c]]
        lambda_avg::Complex{Float64} = Statistics.mean(lambdas)
        index::Int64 = argmin(abs.(lambdas.-lambda_avg))
        v_avg::Vector{Complex{Float64}} = vs[index]
        Lambda_avg[c] = lambda_avg
        V_avg[:,c] = v_avg
    end

    return (Lambda_avg, V_avg)
end

"""
    getBrouckeStability(periodicOrbit)

Return Broucke stability parameters

# Arguments
- `periodicOrbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
"""
function getBrouckeStability(periodicOrbit::CR3BPMSPeriodicOrbit)
    alpha::Float64 = 2-LinearAlgebra.tr(periodicOrbit.monodromy)

    return [alpha, 0.5*((alpha^2)+2-LinearAlgebra.tr(periodicOrbit.monodromy^2))]
end

"""
    getEigenData(periodicOrbit)

Return eigenvalues and -vectors

# Arguments
- `periodicOrbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
"""
function getEigenData(periodicOrbit::CR3BPMSPeriodicOrbit)
    E::LinearAlgebra.Eigen = LinearAlgebra.eigen(periodicOrbit.monodromy)

    return (Vector{Complex{Float64}}(E.values), Matrix{Complex{Float64}}(E.vectors))
end

"""
    getJacobiConstant(periodicOrbit)

Return Jacobi constant

# Arguments
- `periodicOrbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
"""
function getJacobiConstant(periodicOrbit::CR3BPMSPeriodicOrbit)
    return getJacobiConstant(periodicOrbit.dynamicsModel, periodicOrbit.initialCondition)
end

"""
    getStabilityIndex(periodicOrbit)

Return stability index

# Arguments
- `periodicOrbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
"""
function getStabilityIndex(periodicOrbit::CR3BPMSPeriodicOrbit)
    return LinearAlgebra.norm(getEigenData(periodicOrbit)[1], Inf)
end

"""
    getTimeConstant(periodicOrbit)

Return time constant [ndim]

# Arguments
- `periodicOrbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
"""
function getTimeConstant(periodicOrbit::CR3BPMSPeriodicOrbit)
    return periodicOrbit.period/log(getStabilityIndex(periodicOrbit))
end
