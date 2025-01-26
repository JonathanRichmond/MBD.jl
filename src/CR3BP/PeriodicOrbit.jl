"""
CR3BP periodic orbit wrapper

Author: Jonathan Richmond
C: 1/16/23
U: 1/23/25
"""

import DifferentialEquations, LinearAlgebra, StaticArrays
import MBD: CR3BPPeriodicOrbit

export getBrouckeStability, getEigenData, getJacobiConstant, getManifoldArcByTime
export getManifoldByArclength, getManifoldByStepOff, getManifoldByTime, getStabilityIndex
export getTimeConstant

"""
    getBrouckeStability(periodicOrbit)

Return Broucke stability parameters

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getBrouckeStability(periodicOrbit::CR3BPPeriodicOrbit)
    alpha::Float64 = 2-LinearAlgebra.tr(periodicOrbit.monodromy)

    return [alpha, 0.5*((alpha^2)+2-LinearAlgebra.tr(periodicOrbit.monodromy^2))]
end

"""
    getEigenData(periodicOrbit)

Return eigenvalues and -vectors

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getEigenData(periodicOrbit::CR3BPPeriodicOrbit)
    E::LinearAlgebra.Eigen = LinearAlgebra.eigen(periodicOrbit.monodromy)

    return (Vector{Complex{Float64}}(E.values), Matrix{Complex{Float64}}(E.vectors))
end

"""
    getJacobiConstant(periodicOrbit)

Return Jacobi constant

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getJacobiConstant(periodicOrbit::CR3BPPeriodicOrbit)
    return getJacobiConstant(periodicOrbit.dynamicsModel, periodicOrbit.initialCondition)
end

"""
    getManifoldArcByTime(periodicOrbit, stabilitity, direction, d, orbitTime)

Return stable or unstable manifold tubes spaced by time

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `stability::String`: Desired manifold stability
- `direction::String`: Step-off direction
- `d::Float64`: Step-off distance [ndim]
- `orbitTime::Float64`: Normalized time along orbit from initial condition
"""
function getManifoldArcByTime(periodicOrbit::CR3BPPeriodicOrbit, stability::String, direction::String, d::Float64, orbitTime::Float64)
    (eigenvalues::Vector{Complex{Float64}}, eigenvectors::Matrix{Complex{Float64}}) = getEigenData(periodicOrbit)
    index::Int16 = (stability == "Stable") ? Int16(argmin(abs.(eigenvalues))) : Int16(argmax(abs.(eigenvalues)))
    eigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(eigenvectors[:,index])
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    orbitArc::MBD.CR3BPArc = propagate(propagator, appendExtraInitialConditions(periodicOrbit.dynamicsModel, periodicOrbit.initialCondition, MBD.STM), [0, orbitTime*periodicOrbit.period], periodicOrbit.dynamicsModel)
    stateSize::Int64 = getStateSize(periodicOrbit.dynamicsModel, MBD.STM)
    q::StaticArrays.SVector{stateSize, Float64} = StaticArrays.SVector{stateSize, Float64}(getStateByIndex(orbitArc, -1))
    state::Vector{Float64} = q[1:6]
    Phi::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]])
    arcEigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(Phi*eigenvector)
    normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
    step::Int16 = (direction == "Negative") ? Int16(-1) : Int16(1)
    manifoldArc = MBD.CR3BPManifoldArc(periodicOrbit, orbitTime, d, Vector{Complex{Float64}}(state+step*d.*normEigenvector))

    return manifoldArc
end

"""
    getManifoldByArclength(periodicOrbit, stabilitity, direction, d, nArcs)

Return stable or unstable manifold tubes spaced by arclength

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `stability::String`: Desired manifold stability
- `direction::String`: Step-off direction
- `d::Float64`: Step-off distance [ndim]
- `nArcs::Int64`: Number of manifold arcs
"""
function getManifoldByArclength(periodicOrbit::CR3BPPeriodicOrbit, stability::String, direction::String, d::Float64, nArcs::Int64)
    (eigenvalues::Vector{Complex{Float64}}, eigenvectors::Matrix{Complex{Float64}}) = getEigenData(periodicOrbit)
    index::Int16 = (stability == "Stable") ? Int16(argmin(abs.(eigenvalues))) : Int16(argmax(abs.(eigenvalues)))
    eigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(eigenvectors[:,index])
    propagator = MBD.Propagator()
    propagator.equationType = MBD.ARCLENGTH
    orbitArc::MBD.CR3BPArc = propagate(propagator, appendExtraInitialConditions(periodicOrbit.dynamicsModel, periodicOrbit.initialCondition, MBD.ARCLENGTH), [0, periodicOrbit.period], periodicOrbit.dynamicsModel)
    orbitLength::Float64 = getStateByIndex(orbitArc, -1)[43]
    arclength::Vector{Float64} = collect(range(0, orbitLength, nArcs+1))
    manifold = MBD.CR3BPManifold(periodicOrbit, stability, direction)
    stateSize::Int64 = getStateSize(periodicOrbit.dynamicsModel, MBD.ARCLENGTH)
    for a::Int16 in Int16(2):Int16(nArcs+1)
        arclengthEvent = DifferentialEquations.ContinuousCallback(arclengthCondition, terminateAffect!)
        arc::MBD.CR3BPArc = propagateWithEvent(propagator, arclengthEvent, appendExtraInitialConditions(periodicOrbit.dynamicsModel, periodicOrbit.initialCondition, MBD.ARCLENGTH), [0, periodicOrbit.period], periodicOrbit.dynamicsModel, [arclength[a]])
        q::StaticArrays.SVector{stateSize, Float64} = StaticArrays.SVector{stateSize, Float64}(getStateByIndex(arc, -1))
        t::Float64 = getTimeByIndex(arc, -1)
        state::Vector{Float64} = q[1:6]
        Phi::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]])
        arcEigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(Phi*eigenvector)
        normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
        step::Int16 = (direction == "Negative") ? Int16(-1) : Int16(1)
        push!(manifold.initialConditions, state+step*d.*normEigenvector)
        push!(manifold.orbitTimes, t/periodicOrbit.period)
        push!(manifold.ds, d)
    end

    return manifold
end

"""
    getManifoldByStepOff(periodicOrbit, stability, direction, d_max, nArcs)

Return stable or unstable manifold tubes spaced by step-off

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `stability::String`: Desired manifold stability
- `direction::String`: Step-off direction
- `d_max::Float64`: Maximum step-off distance [ndim]
- `nArcs::int64`: Number of manifold arcs
"""
function getManifoldByStepOff(periodicOrbit::CR3BPPeriodicOrbit, stability::String, direction::String, d_max::Float64, nArcs::Int64)
    (eigenvalues::Vector{Complex{Float64}}, eigenvectors::Matrix{Complex{Float64}}) = getEigenData(periodicOrbit)
    index::Int16 = (stability == "Stable") ? Int16(argmin(abs.(eigenvalues))) : Int16(argmax(abs.(eigenvalues)))
    eigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(eigenvectors[:,index])
    normEigenvector::Vector{Complex{Float64}} = eigenvector./LinearAlgebra.norm(eigenvector[1:3])
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    lambda::Float64 = eigenvalues[index]
    manifold = MBD.CR3BPManifold(periodicOrbit, stability, direction)
    for a::Int16 in Int16(1):Int16(nArcs)
        beta::Float64 = (direction == "Negative") ? (pi-acos(lambda^((a-1)/(nArcs-1)-1))) : acos(lambda^((a-1)/(nArcs-1)-1))
        d::Float64 = d_max*cos(beta)
        push!(manifold.initialConditions, periodicOrbit.initialCondition+d.*normEigenvector)
        push!(manifold.orbitTimes, 0.0)
        push!(manifold.ds, d)
    end

    return manifold
end

"""
    getManifoldByTime(periodicOrbit, stabilitity, direction, d, nArcs)

Return stable or unstable manifold tubes spaced by time

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `stability::String`: Desired manifold stability
- `direction::String`: Step-off direction
- `d::Float64`: Step-off distance [ndim]
- `nArcs::Int64`: Number of manifold arcs
"""
function getManifoldByTime(periodicOrbit::CR3BPPeriodicOrbit, stability::String, direction::String, d::Float64, nArcs::Int64)
    (eigenvalues::Vector{Complex{Float64}}, eigenvectors::Matrix{Complex{Float64}}) = getEigenData(periodicOrbit)
    index::Int16 = (stability == "Stable") ? Int16(argmin(abs.(eigenvalues))) : Int16(argmax(abs.(eigenvalues)))
    eigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(eigenvectors[:,index])
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    time::Vector{Float64} = collect(range(0, periodicOrbit.period, nArcs+1))
    manifold = MBD.CR3BPManifold(periodicOrbit, stability, direction)
    stateSize::Int64 = getStateSize(periodicOrbit.dynamicsModel, MBD.STM)
    for a::Int16 in Int16(2):Int16(nArcs+1)
        arc::MBD.CR3BPArc = propagate(propagator, appendExtraInitialConditions(periodicOrbit.dynamicsModel, periodicOrbit.initialCondition, MBD.STM), [0, time[a]], periodicOrbit.dynamicsModel)
        q::StaticArrays.SVector{stateSize, Float64} = StaticArrays.SVector{stateSize, Float64}(getStateByIndex(arc, -1))
        state::Vector{Float64} = q[1:6]
        Phi::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([q[7:12] q[13:18] q[19:24] q[25:30] q[31:36] q[37:42]])
        arcEigenvector::StaticArrays.SVector{6, Complex{Float64}} = StaticArrays.SVector{6, Complex{Float64}}(Phi*eigenvector)
        normEigenvector::Vector{Complex{Float64}} = arcEigenvector./LinearAlgebra.norm(arcEigenvector[1:3])
        step::Int16 = (direction == "Negative") ? Int16(-1) : Int16(1)
        push!(manifold.initialConditions, state+step*d.*normEigenvector)
        push!(manifold.orbitTimes, time[a]/periodicOrbit.period)
        push!(manifold.ds, d)
    end

    return manifold
end

"""
    getStabilityIndex(periodicOrbit)

Return stability index

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getStabilityIndex(periodicOrbit::CR3BPPeriodicOrbit)
    return LinearAlgebra.norm(getEigenData(periodicOrbit)[1], Inf)
end

"""
    getTimeConstant(periodicOrbit)

Return time constant [ndim]

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getTimeConstant(periodicOrbit::CR3BPPeriodicOrbit)
    return periodicOrbit.period/log(getStabilityIndex(periodicOrbit))
end
