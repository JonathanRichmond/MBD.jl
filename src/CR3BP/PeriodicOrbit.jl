"""
CR3BP periodic orbit wrapper

Author: Jonathan Richmond
C: 1/16/23
U: 9/10/23
"""

import DifferentialEquations, LinearAlgebra
import MBD: CR3BPPeriodicOrbit

export getManifold, getStability!

"""
    getManifold(periodicOrbit, dynamicsModel, stabilitity, d)

Return stable or unstable manifold tubes

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model
- `stability::String`: Desired manifold stability
- `d::Float64`: Stepoff distance [ndim]
"""
function getManifold(periodicOrbit::CR3BPPeriodicOrbit, dynamicsModel::CR3BPDynamicsModel, stability::String, d::Float64)
    index::Int64 = (stability == "Stable") ?  argmax(abs.(periodicOrbit.eigenvalues)) : argmin(abs.(periodicOrbit.eigenvalues))
    eigenvector::Vector{Complex{Float64}} = periodicOrbit.eigenvectors[:,index]
    propagator = MBD.Propagator()
    propagator.equationType = MBD.ARCLENGTH
    orbitArc::MBD.Arc = propagate(propagator, vcat(appendExtraInitialConditions(dynamicsModel, periodicOrbit.initialCondition, MBD.STM), 0.0), [0, periodicOrbit.period], dynamicsModel)
    orbitLength::Float64 = getStateByIndex(orbitArc, -1)[43]
    arclength::Vector{Float64} = collect(range(0, orbitLength, 11))
    posManifold::Vector{MBD.CR3BPManifoldArc} = []
    negManifold::Vector{MBD.CR3BPManifoldArc} = []
    for a::Int64 in 2:length(arclength)
        callbackEvent = DifferentialEquations.ContinuousCallback(arclengthCondition, terminateAffect!)
        arc::MBD.Arc = propagateWithEvent(propagator, callbackEvent, vcat(appendExtraInitialConditions(dynamicsModel, periodicOrbit.initialCondition, MBD.STM), 0.0), [0, periodicOrbit.period], dynamicsModel, [arclength[a]])
        q::Vector{Float64} = getStateByIndex(arc, -1)
        state::Vector{Float64} = q[1:6]
        Phi::Matrix{Float64} = [q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]]
        arcEigenvector::Vector{Complex{Float64}} = Phi*eigenvector
        normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
        posManifoldArc = MBD.CR3BPManifoldArc(state+d.*normEigenvector, periodicOrbit)
        negManifoldArc = MBD.CR3BPManifoldArc(periodicOrbit.initialCondition-d.*normEigenvector, periodicOrbit)
        push!(posManifold, posManifoldArc)
        push!(negManifold, negManifoldArc)
    end

    return (posManifold, negManifold)
end

"""
    getStability!(periodicOrbit)

Return periodic orbit object with updated stability properties

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getStability!(periodicOrbit::CR3BPPeriodicOrbit)
    E::LinearAlgebra.Eigen = LinearAlgebra.eigen(periodicOrbit.monodromy)
    periodicOrbit.eigenvalues = E.values
    periodicOrbit.eigenvectors = E.vectors
    alpha::Float64 = 2-LinearAlgebra.tr(periodicOrbit.monodromy)
    beta::Float64 = 0.5*((alpha^2)+2-LinearAlgebra.tr(periodicOrbit.monodromy^2))
    periodicOrbit.BrouckeStability = [alpha, beta]
    periodicOrbit.nu = LinearAlgebra.norm(periodicOrbit.eigenvalues, Inf)
    periodicOrbit.tau = periodicOrbit.period/log(periodicOrbit.nu)
end

function arclengthCondition(state, time, integrator)
    state[43]-integrator.p[2][1]
end

function terminateAffect!(integrator)
    DifferentialEquations.terminate!(integrator)
end
