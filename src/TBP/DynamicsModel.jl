"""
TBP dynamics model wrapper

Author: Jonathan Richmond
C: 9/14/23
"""

import LinearAlgebra
import MBD: TBPDynamicsModel

export appendExtraInitialConditions, evaluateEquations, getEquationsOfMotion
export getEpochDependencies, getLambertArc, getParameterDependencies
export getPrimaryPosition, getStateSize, getStateTransitionMatrix
export isEpochIndependent

"""
    appendExtraInitialConditions(dynamicsModel, q0_simple, outputEquationType)

Return state vector with extra initial conditions

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `q0_simple::Vector{Float64}`: Simple initial state vector [ndim]
- `outputEquationType::EquationType`: Output state EOM type
"""
function appendExtraInitialConditions(dynamicsModel::TBPDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::MBD.EquationType)
    n_in::Int64 = getStateSize(dynamicsModel, MBD.SIMPLE)
    (length(q0_simple) == n_in) || throw(ArgumentError("State vector length is $(length(q0_simple)), but should be $n_in"))
    n_simple::Int64 = getStateSize(dynamicsModel, MBD.SIMPLE)
    n_STM::Int64 = getStateSize(dynamicsModel, MBD.STM)
    n_out::Int64 = getStateSize(dynamicsModel, outputEquationType)
    q0_out::Vector{Float64} = zeros(Float64, n_out)
    if n_in >= n_out
        q0_out = q0_simple[1:n_out]
    else
        q0_out[1:n_in] = q0_simple
        (n_in == n_simple) && [q0_out[i] = 1 for i in n_simple+1:n_simple+1:n_STM]
    end

    return q0_out
end

"""
    evaluateEquations(dynamicsModel, equationType, t, q)

Return time derivative of state vector

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `equationType::EquationType`: EOM type
- `t::Float64`: Time [ndim]
- `q::Vector{Float64}`: State vector [ndim]
"""
function evaluateEquations(dynamicsModel::TBPDynamicsModel, equationType::MBD.EquationType, t::Float64, q::Vector{Float64})
    qdot::Vector{Float64} = Vector{Float64}(undef, getStateSize(dynamicsModel, equationType))
    EOMs::MBD.TBPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, equationType)
    computeDerivatives!(qdot, q, EOMs, t)

    return qdot
end

"""
    getEquationsOfMotion(dynamicsModel, equationType)

Return EOMs

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `equationType::EquationType`: EOM type
"""
function getEquationsOfMotion(dynamicsModel::TBPDynamicsModel, equationType::MBD.EquationType)
    return MBD.TBPEquationsOfMotion(equationType, dynamicsModel)
end

"""
    getEpochDependencies(dynamicsModel, q)

Return derivative of state with respect to epoch

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getEpochDependencies(dynamicsModel::TBPDynamicsModel, q_full::Vector{Float64})
    n_full::Int64 = getStateSize(dynamicsModel, MBD.FULL)
    (length(q_full) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int64 = getStateSize(dynamicsModel, MBD.SIMPLE)

    isEpochIndependent(dynamicsModel) ? (return zeros(Float64, n_simple)) : (return q_full[n_simple*(n_simple+1)+1:n_simple*(n_simple+1)+n_simple])
end

"""
    getLambertArc(initialPos, finalPos, TOF, transferMethod)

Return initial/final velocities

# Arguments
- `initalPos::Vector{Float64}`: Initial position [ndim]
- `finalPos::Vector{Float64}`: Final position [ndim]
- `TOF::Float64`: Time of flight [ndim]
- `transferMethod::String`: Transfer method
"""
function getLambertArc(initialPos::Vector{Float64}, finalPos::Vector{Float64}, TOF::Float64, transferMethod::String)
    r0::Float64 = LinearAlgebra.norm(initialPos)
    rf::Float64 = LinearAlgebra.norm(finalPos)
    cosdeltanu::Float64 = LinearAlgebra.dot(initialPos, finalPos)/(r0*rf)
    t_m::Float64 = (transferMethod == "Short") ? 1.0 : -1.0
    A::Float64 = t_m*sqrt(rf*r0*(1+cosdeltanu))
    (A == 0) && throw(ErrorException("A = 0 so Lamber arc cannot be computed"))
    psi_n::Float64 = 0.0
    c_2::Float64 = 1/2
    c_3::Float64 = 1/6
    psi_up::Float64 = 4*pi^2
    psi_low::Float64 = -4*pi
    deltat_n::Float64 = TOF+100
    while abs(deltat_n-TOF) >= 1E-9
        y_n::Float64 = r0+rf+(A*(psi_n*c_3-1))/sqrt(c_2)
        x_n::Float64 = sqrt(y_n/c_2)
        deltat_n = x_n^3*c_3+A*sqrt(y_n)
        (deltat_n > deltat) ? (psi_up = psi_n) : (psi_low = psi_n)
        psi_n = (psi_up+psi_low)/2
        if psi_n > 1E-6
            c_2 = (1-cos(sqrt(psi_n)))/psi_n
            c_3 = (sqrt(psi_n)-sin(sqrt(psi_n)))/sqrt(psi_n^3)
        elseif psi_n < -1E-6
            c_2 = (1-cosh(sqrt(-psi_n)))/psi_n
            c_3 = (sinh(sqrt(-psi_n))-sqrt(-psi_n))/sqrt((-psi_n)^3)
        else
            c_2 = 1/2
            c_3 = 1/6
        end
    end
    f::Float64 = 1-y_n/r0
    gdot::Float64 = 1-y_n/rf
    g::Float64 = A*sqrt(y_n)

    return ((finalPos-f.*initialPos)./g, (gdot.*finalPos-initialPos)./g)
end

"""
    getParameterDependencies(dynamicsModel, q_full)

Return derivative of state with respect to parameters

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getParameterDependencies(dynamicsModel::TBPDynamicsModel, q_full::Vector{Float64})
    n_full::Int64 = getStateSize(dynamicsModel, MBD.FULL)
    (length(q_full) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int64 = getStateSize(dynamicsModel, MBD.SIMPLE)
    i0::Int64 = n_simple*(n_simple+1+(isEpochIndependent(dynamicsModel) ? 0 : 1))
    n::Int64 = n_full-i0
    (n == 0) && (return zeros(Float64, (n_simple, 0)))
    n_params::Int64 = n/n_simple
    dqdp::Matrix{Float64} = zeros(Float64, (n_simple, n_params))
    for r::Int64 in 1:n_simple
        for c::Int64 in 1:n_params
            dqdp[r,c] = q_full[i0+n_simple*(c-1)+r]
        end
    end

    return dqdp
end

"""
    getPrimaryPosition(dynamicsModel)

Return location of primary

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
"""
function getPrimaryPosition(dynamicsModel::TBPDynamicsModel)
    return [0.0, 0.0, 0.0]
end

"""
    getStateSize(dynamicsModel, equationType)

Return number of state variables

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `equationType::EquationType`: EOM type
"""
function getStateSize(dynamicsModel::TBPDynamicsModel, equationType::MBD.EquationType)
    type = Dict(MBD.SIMPLE => 6, MBD.STM => 42, MBD.FULL => 42)

    return type[equationType]
end

"""
    getStateTransitionMatrix(dynamicsModel, q0_STM)

Return STM

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `q0_STM::Vector{Float64}`: Initial state vector with STM in column-major order [ndim]
"""
function getStateTransitionMatrix(dynamicsModel::TBPDynamicsModel, q0_STM::Vector{Float64})
    n_STM::Int64 = getStateSize(dynamicsModel, MBD.STM)
    (length(q0_STM) < n_STM) && throw(ArgumentError("State vector length is $(length(q0_STM)), but should be $n_STM"))
    n_simple::Int64 = getStateSize(dynamicsModel, MBD.SIMPLE)
    STM::Matrix{Float64} = zeros(Float64, (n_simple, n_simple))
    for r::Int64 in 1:n_simple
        for c::Int64 in 1:n_simple
            STM[r,c] = q0_STM[n_simple+n_simple*(c-1)+r]
        end
    end

    return STM
end

"""
    isEpochIndependent(dynamicsModel)

Return true if dynamics model is epoch independent

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
"""
function isEpochIndependent(dynamicsModel::TBPDynamicsModel)
    return true
end
