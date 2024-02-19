"""
TBP dynamics model wrapper

Author: Jonathan Richmond
C: 9/14/23
"""

import LinearAlgebra
import MBD: TBPDynamicsModel

export appendExtraInitialConditions, evaluateEquations, getCartesianState
export getEquationsOfMotion, getEpochDependencies, getExcursion, getLambertArc
export getOsculatingOrbitalElements, getParameterDependencies, getPeriod
export getPrimaryPosition, getResonantOrbit, getStateSize
export getStateTransitionMatrix, isEpochIndependent, primaryInertial2Rotating
export solveKeplersEquation

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
    computeDerivatives!(qdot, q, (EOMs,), t)

    return qdot
end

"""
    getCartesianState(dynamicsModel, a, e, i, Omega, omega, theta)

Return 2BP trajectory object with Cartesian state

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `a::Float64`: Semimajor axis [km]
- `e::Float64`: Eccentricity
- `i::Float64`: Inclination [rad]
- `Omega::Float64`: Longitude of ascending node [rad]
- `omega::Float64`: Argument of periapsis [rad]
- `theta::Float64`: True anomaly [rad]
"""
function getCartesianState(dynamicsModel::TBPDynamicsModel, a::Float64, e::Float64, i::Float64, Omega::Float64, omega::Float64, theta::Float64)
    eccentricAnomaly::Float64 = atan(sqrt(1-e^2)*sin(theta), e+cos(theta))
    circularRadius::Float64 = a*(1-e*cos(eccentricAnomaly))
    r_E::Vector{Float64} = [circularRadius*cos(theta), circularRadius*sin(theta), 0]
    v_E::Vector{Float64} = (sqrt(dynamicsModel.systemData.gravParam*a)/circularRadius).*[-sin(eccentricAnomaly), sqrt(1-e^2)*cos(eccentricAnomaly), 0]
    C::Matrix{Float64} = [cos(omega)*cos(Omega)-sin(omega)*cos(i)*sin(Omega) -sin(omega)*cos(Omega)-cos(omega)*cos(i)*sin(Omega) 0; cos(omega)*sin(Omega)+sin(omega)*cos(i)*cos(Omega) -sin(omega)*sin(Omega)+cos(omega)*cos(i)*cos(Omega) 0; sin(omega)*sin(i) cos(omega)*sin(i) 0]
    pos_dim::Vector{Float64} = C*r_E
    vel_dim::Vector{Float64} = C*v_E
    stateTrajectory = MBD.TBPTrajectory(append!(copy(pos_dim), vel_dim), dynamicsModel)
    angularMomentum::Vector{Float64} = LinearAlgebra.cross(pos_dim, vel_dim)
    stateTrajectory.h = LinearAlgebra.norm(angularMomentum)
    stateTrajectory.i = i
    stateTrajectory.Omega = Omega
    stateTrajectory.e = e
    stateTrajectory.a = a
    stateTrajectory.r_p = a*(1-e)
    stateTrajectory.r_a = a*(1+e)
    stateTrajectory.omega = omega
    stateTrajectory.theta = theta
    stateTrajectory.E = -dynamicsModel.systemData.gravParam/(2*a)

    return stateTrajectory
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
    getExcursion(dynamicsModel, q)

Return distance from primary

# Arguments
- `dynamicsModel::TBPDynamicsModel`: CR3BP dynamics model object
- `q::Vector{Float64}`: State vector [ndim]
"""
function getExcursion(dynamicsModel::TBPDynamicsModel, q::Vector{Float64})
    primaryPos::Vector{Float64} = getPrimaryPosition(dynamicsModel)

    return LinearAlgebra.norm(q[1:3]-primaryPos)*dynamicsModel.systemData.charLength
end

"""
    getLambertArc(dynamicsModel, initialPos, finalPos, TOF, transferMethod)

Return initial/final velocities

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `initalPos::Vector{Float64}`: Initial position [dim]
- `finalPos::Vector{Float64}`: Final position [dim]
- `TOF::Float64`: Time of flight [dim]
- `transferMethod::String`: Transfer method
"""
function getLambertArc(dynamicsModel::TBPDynamicsModel, initialPos::Vector{Float64}, finalPos::Vector{Float64}, TOF::Float64, transferMethod::String)
    r0::Float64 = LinearAlgebra.norm(initialPos)
    rf::Float64 = LinearAlgebra.norm(finalPos)
    cosdeltanu::Float64 = LinearAlgebra.dot(initialPos, finalPos)/(r0*rf)
    t_m::Float64 = (transferMethod == "Short") ? 1.0 : -1.0
    A::Float64 = t_m*sqrt(rf*r0*(1+cosdeltanu))
    (A == 0) && throw(ErrorException("A = 0 so Lambert arc cannot be computed"))
    psi_n::Float64 = 0.0
    c_2::Float64 = 1/2
    c_3::Float64 = 1/6
    psi_up::Float64 = 4*pi^2
    psi_low::Float64 = -8*pi
    deltat_n::Float64 = TOF+100
    iter::Int64 = 0
    while (abs(deltat_n-TOF) >= 1E-6) && (iter < 200)
        y_n::Float64 = r0+rf+(A*(psi_n*c_3-1))/sqrt(c_2)
        if (A > 0) && (y_n < 0)
            psi_low /= 2
        else
            x_n::Float64 = sqrt(y_n/c_2)
            deltat_n = (x_n^3*c_3+A*sqrt(y_n))/sqrt(dynamicsModel.systemData.gravParam)
            (deltat_n > TOF) ? (psi_up = psi_n) : (psi_low = psi_n)
        end
        iter += 1
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
    (iter < 100) || throw(ErrorException("Could not converge Lambert arc"))
    y_n = r0+rf+(A*(psi_n*c_3-1))/sqrt(c_2)
    f::Float64 = 1-y_n/r0
    gdot::Float64 = 1-y_n/rf
    g::Float64 = A*sqrt(y_n/dynamicsModel.systemData.gravParam)

    return ((finalPos-f.*initialPos)./g, (gdot.*finalPos-initialPos)./g)
end

"""
    getOsculatingOrbitalElements(dynamicsModel, state_dim)

Return 2BP trajectory object with orbital elements

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `state_dim::Vector{Float64}`: Primary-centered inertial state [dim]
"""
function getOsculatingOrbitalElements(dynamicsModel::TBPDynamicsModel, state_dim::Vector{Float64})
    trajectory = MBD.TBPTrajectory(state_dim, dynamicsModel)
    r::Float64 = LinearAlgebra.norm(state_dim[1:3])
    v_r::Float64 = LinearAlgebra.dot(state_dim[4:6], state_dim[1:3])/r
    angularMomentum::Vector{Float64} = LinearAlgebra.cross(state_dim[1:3], state_dim[4:6])
    trajectory.h = LinearAlgebra.norm(angularMomentum)
    trajectory.i = acos(angularMomentum[3]/trajectory.h)
    n::Vector{Float64} = LinearAlgebra.cross([0, 0, 1], angularMomentum)
    trajectory.Omega = (n[2] < 0) ? 2*pi-acos(n[1]/LinearAlgebra.norm(n)) : acos(n[1]/LinearAlgebra.norm(n))
    eccentricityVector::Vector{Float64} = LinearAlgebra.cross(state_dim[4:6], angularMomentum)./dynamicsModel.systemData.gravParam-state_dim[1:3]./r
    trajectory.e = LinearAlgebra.norm(eccentricityVector)
    trajectory.a = trajectory.h^2/(dynamicsModel.systemData.gravParam*(1-trajectory.e^2))
    trajectory.r_p = trajectory.a*(1-trajectory.e)
    trajectory.r_a = trajectory.a*(1+trajectory.e)
    trajectory.omega = (eccentricityVector[3] < 0) ? 2*pi-acos(LinearAlgebra.dot(n, eccentricityVector)/(LinearAlgebra.norm(n)*trajectory.e)) : acos(LinearAlgebra.dot(n, eccentricityVector)/(LinearAlgebra.norm(n)*trajectory.e))
    trajectory.theta = (v_r < 0) ? 2*pi-acos(LinearAlgebra.dot(eccentricityVector, state_dim[1:3])/(trajectory.e*r)) : acos(LinearAlgebra.dot(eccentricityVector, state_dim[1:3])/(trajectory.e*r))
    trajectory.E = -dynamicsModel.systemData.gravParam/(2*trajectory.a)

    return trajectory
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
    getPeriod(dynamicsModel, trajectory)

Return two-body orbital period

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `trajectory::TBPTrajectory`: TBP trajectory object
"""
function getPeriod(dynamicsModel::TBPDynamicsModel, trajectory::MBD.TBPTrajectory)
    meanMotion::Float64 = sqrt(dynamicsModel.systemData.gravParam/abs(trajectory.a)^3)

    return 2*pi/meanMotion
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
    getResonantOrbit(dynamicsModel, secondaryData, p, q, e)

Return two-body resonant orbit and period

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `secondaryData::BodyData`: Body data object
- `p::Int64`: Spacecraft revolutions around primary
- `q::Int64`: Secondary revolutions around primary
- `e::Float64`: Eccentricity
"""
function getResonantOrbit(dynamicsModel::TBPDynamicsModel, secondaryData::MBD.BodyData, p::Int64, q::Int64, e::Float64)
    a_p::Float64 = (q^2*secondaryData.orbitRadius^3/(p^2))^(1/3)
    r_p_p::Float64 = a_p*(1-e)
    v::Float64 = sqrt(2*dynamicsModel.systemData.gravParam*(1/r_p_p-1/(2*a_p)))
    q0Dim::Vector{Float64} = [r_p_p, 0, 0, 0, v, 0]
    resonantOrbit = MBD.TBPTrajectory(q0Dim, dynamicsModel)
    resonantOrbit.a = a_p
    resonantOrbit.e = e

    return (resonantOrbit, getPeriod(dynamicsModel, resonantOrbit)*p)
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

"""
    primaryInertial2Rotating(dynamicsModel, secondaryData, states_inertial, times)

Return rotating frame states

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `secondaryData::BodyData`: Body data object
- `states_primaryInertial::Vector{Vector{Float64}}`: Primary-centered inertial states [dim]
- `times::Vector{Float64}`: Epochs [dim]
"""
function primaryInertial2Rotating(dynamicsModel::TBPDynamicsModel, secondaryData::MBD.BodyData, states_primaryInertial::Vector{Vector{Float64}}, times::Vector{Float64})
    (length(states_primaryInertial) == length(times)) || throw(ArgumentError("Number of state vectors, $(length(states_primaryInertial)), must match number of times, $(length(times))"))
    T::Float64 = 2*pi*secondaryData.orbitRadius^(3/2)/sqrt(dynamicsModel.systemData.gravParam)
    t::Vector{Float64} = 2*pi*times./T
    states::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
    for i in 1:length(times)
        C::Matrix{Float64} = [cos(t[i]) -sin(t[i]) 0; sin(t[i]) cos(t[i]) 0; 0 0 1]
        Cdot::Matrix{Float64} = (2*pi/T)*[-sin(t[i]) -cos(t[i]) 0; cos(t[i]) -sin(t[i]) 0; 0 0 0]
        N::Matrix{Float64} = [C zeros(Float64, (3,3)); Cdot C]
        states[i] = N\states_primaryInertial[i]
    end

    return states
end

"""
    solveKeplersEquation(dynamicsModel, trajectory)

Return time since periapsis based on true anomaly

# Arguments
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
- `trajectory::TBPTrajectory`: TBP trajectory object
"""
function solveKeplersEquation(dynamicsModel::TBPDynamicsModel, trajectory::MBD.TBPTrajectory)
    meanMotion::Float64 = sqrt(dynamicsModel.systemData.gravParam/abs(trajectory.a)^3)
    (trajectory.theta == 1.0*pi) && (return pi/meanMotion)
    eccentricAnomaly::Float64 = 2*atan(tan(trajectory.theta/2)/sqrt((1+trajectory.e)/(1-trajectory.e)))
    meanAnomaly::Float64 = eccentricAnomaly-trajectory.e*sin(eccentricAnomaly)
    meanAnomaly < 0 ? (return (meanAnomaly+2*pi)/meanMotion) : (return meanAnomaly/meanMotion)
end
