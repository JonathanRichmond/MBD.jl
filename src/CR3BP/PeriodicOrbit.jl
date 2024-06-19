"""
CR3BP periodic orbit wrapper

Author: Jonathan Richmond
C: 1/16/23
U: 9/10/23
"""

import DifferentialEquations, LinearAlgebra
import MBD: CR3BPPeriodicOrbit

export getManifoldArcByTime, getManifoldByArclength, getManifoldByTime
export getStability!

"""
    getManifoldArcByTime(periodicOrbit, dynamicsModel, stabilitity, d, orbitTime, direction)

Return stable or unstable manifold tubes spaced by time

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model
- `stability::String`: Desired manifold stability
- `d::Float64`: Step-off distance [ndim]
- `orbitTime::Float64`: Time along orbit from initial condition [ndim]
- `direction::String`: Step-off direction
"""
function getManifoldArcByTime(periodicOrbit::CR3BPPeriodicOrbit, dynamicsModel::CR3BPDynamicsModel, stability::String, d::Float64, orbitTime::Float64, directino::String)
    index::Int64 = (stability == "Stable") ?  argmin(abs.(periodicOrbit.eigenvalues)) : argmax(abs.(periodicOrbit.eigenvalues))
    eigenvector::Vector{Complex{Float64}} = periodicOrbit.eigenvectors[:,index]
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    orbitArc::MBD.Arc = propagate(propagator, appendExtraInitialConditions(dynamicsModel, periodicOrbit.initialCondition, MBD.STM), [0, orbitTime], dynamicsModel)
    q::Vector{Float64} = getStateByIndex(orbitArc, -1)
    state::Vector{Float64} = q[1:6]
    Phi::Matrix{Float64} = [q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]]
    arcEigenvector::Vector{Complex{Float64}} = Phi*eigenvector
    normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
    step::Int64 = (direction == "Negative") ? -1 : 1
    manifoldArc = MBD.CR3BPManifoldArc(state+step*d.*normEigenvector, periodicOrbit, orbitTime)

    return manifoldArc
end

"""
    getManifoldByArclength(periodicOrbit, dynamicsModel, stabilitity, d, nArcs)

Return stable or unstable manifold tubes spaced by arclength

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model
- `stability::String`: Desired manifold stability
- `d::Float64`: Step-off distance [ndim]
- `nArcs::Int64`: Number of manifold arcs
"""
function getManifoldByArclength(periodicOrbit::CR3BPPeriodicOrbit, dynamicsModel::CR3BPDynamicsModel, stability::String, d::Float64, nArcs::Int64)
    index::Int64 = (stability == "Stable") ?  argmin(abs.(periodicOrbit.eigenvalues)) : argmax(abs.(periodicOrbit.eigenvalues))
    eigenvector::Vector{Complex{Float64}} = periodicOrbit.eigenvectors[:,index]
    propagator = MBD.Propagator()
    propagator.equationType = MBD.ARCLENGTH
    orbitArc::MBD.Arc = propagate(propagator, vcat(appendExtraInitialConditions(dynamicsModel, periodicOrbit.initialCondition, MBD.STM), 0.0), [0, periodicOrbit.period], dynamicsModel)
    orbitLength::Float64 = getStateByIndex(orbitArc, -1)[43]
    arclength::Vector{Float64} = collect(range(0, orbitLength, nArcs+1))
    posManifold::Vector{MBD.CR3BPManifoldArc} = []
    negManifold::Vector{MBD.CR3BPManifoldArc} = []
    for a::Int64 in 2:nArcs+1
        arclengthEvent = DifferentialEquations.ContinuousCallback(arclengthCondition, terminateAffect!)
        arc::MBD.Arc = propagateWithEvent(propagator, arclengthEvent, vcat(appendExtraInitialConditions(dynamicsModel, periodicOrbit.initialCondition, MBD.STM), 0.0), [0, periodicOrbit.period], dynamicsModel, [arclength[a]])
        q::Vector{Float64} = getStateByIndex(arc, -1)
        t::Float64 = getTimeByIndex(arc, -1)
        state::Vector{Float64} = q[1:6]
        Phi::Matrix{Float64} = [q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]]
        arcEigenvector::Vector{Complex{Float64}} = Phi*eigenvector
        normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
        posManifoldArc = MBD.CR3BPManifoldArc(state+d.*normEigenvector, periodicOrbit, t)
        negManifoldArc = MBD.CR3BPManifoldArc(state-d.*normEigenvector, periodicOrbit, t)
        push!(posManifold, posManifoldArc)
        push!(negManifold, negManifoldArc)
    end

    return (posManifold, negManifold)
end

"""
    getManifoldByTime(periodicOrbit, dynamicsModel, stabilitity, d, nArcs)

Return stable or unstable manifold tubes spaced by time

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model
- `stability::String`: Desired manifold stability
- `d::Float64`: Step-off distance [ndim]
- `nArcs::Int64`: Number of manifold arcs
"""
function getManifoldByTime(periodicOrbit::CR3BPPeriodicOrbit, dynamicsModel::CR3BPDynamicsModel, stability::String, d::Float64, nArcs::Int64)
    index::Int64 = (stability == "Stable") ?  argmin(abs.(periodicOrbit.eigenvalues)) : argmax(abs.(periodicOrbit.eigenvalues))
    eigenvector::Vector{Complex{Float64}} = periodicOrbit.eigenvectors[:,index]
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    time::Vector{Float64} = collect(range(0, periodicOrbit.period, nArcs+1))
    posManifold::Vector{MBD.CR3BPManifoldArc} = []
    negManifold::Vector{MBD.CR3BPManifoldArc} = []
    for a::Int64 in 2:nArcs+1
        arc::MBD.Arc = propagate(propagator, appendExtraInitialConditions(dynamicsModel, periodicOrbit.initialCondition, MBD.STM), [0, time[a]], dynamicsModel)
        q::Vector{Float64} = getStateByIndex(arc, -1)
        state::Vector{Float64} = q[1:6]
        Phi::Matrix{Float64} = [q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]]
        arcEigenvector::Vector{Complex{Float64}} = Phi*eigenvector
        normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
        posManifoldArc = MBD.CR3BPManifoldArc(state+d.*normEigenvector, periodicOrbit, time[a])
        negManifoldArc = MBD.CR3BPManifoldArc(state-d.*normEigenvector, periodicOrbit, time[a])
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

"""
    shallowClone(periodicOrbit)

Return copy of periodic orbit object

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function shallowClone(periodicOrbit::CR3BPPeriodicOrbit)
    object = CR3BPPeriodicOrbit(deepClone(periodicOrbit.problem), periodicOrbit.targeter)
    object.BrouckeStability = copy(periodicOrbit.BrouckeStability)
    object.eigenvalues = copy(periodicOrbit.eigenvalues)
    object.eigenvectors = copy(periodicOrbit.eigenvectors)
    object.initialCondition = copy(periodicOrbit.initialCondition)
    object.JacobiConstant = periodicOrbit.JacobiConstant
    object.monodromy = copy(periodicOrbit.monodromy)
    object.nu = periodicOrbit.nu
    object.tau = periodicOrbit.tau

    return object
end
