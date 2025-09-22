"""
Keplerian dynamics model wrapper

Author: Jonathan Richmond
C: 9/18/25
"""

import LinearAlgebra, SPICE, StaticArrays
import MBD: KDynamicsModel

export appendExtraInitialConditions, checkSTM, evaluateEquations, getCartesianState, getEnergy
export getEpochDependencies, getEquationsOfMotion, getExcursion
export getMassParameter, getOrbitalElements, getParameterDependencies, getPeriod
export getPrimaryState, getStateSize, getStateTransitionMatrix
export isEpochIndependent, solveKeplersEquation

"""
    appendExtraInitialConditions(dynamicsModel, q0_simple, outputEquationType)

Return state vector with extra initial conditions

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `q0_simple::Vector{Float64}`: Simple initial state vector [ndim]
- `outputEquationType::EquationType`: Output state EOM type
"""
function appendExtraInitialConditions(dynamicsModel::KDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::MBD.EquationType)
    n_in::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    (Int16(length(q0_simple)) == n_in) || throw(ArgumentError("State vector length is $(length(q0_simple)), but should be $n_in"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    n_STM::Int16 = getStateSize(dynamicsModel, MBD.STM)
    n_out::Int16 = getStateSize(dynamicsModel, outputEquationType)
    q0_out::Vector{Float64} = zeros(Float64, n_out)
    if n_in >= n_out
        q0_out = q0_simple[1:n_out]
    else
        q0_out[1:n_in] = q0_simple
        [q0_out[i] = 1 for i in n_simple+1:n_simple+1:n_STM]
    end

    return q0_out
end

"""
    checkSTM(dynamicsModel; relTol)

Return true if STM is accurate

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `relTol::Float64`: Relative tolerance (default = 2E-3)
"""
function checkSTM(dynamicsModel::KDynamicsModel, relTol::Float64 = 2E-3)
    stepSize::Float64 = sqrt(eps(Float64))
    numStates::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    propagator = MBD.Propagator()
    propagatorSTM = MBD.Propagator(equationType = MBD.STM)
    X::Vector{Float64} = [400000.0, 0, 0, 0, 1.0, 0]
    tau::Float64 = 3600
    arc::MBD.KArc = propagate(propagatorSTM, appendExtraInitialConditions(dynamicsModel, X, MBD.STM), [0, tau], dynamicsModel)
    STMAnalytical::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64}(getStateTransitionMatrix(dynamicsModel, getStateByIndex(arc, -1)))
    STMNumerical::StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64}(zeros(Float64, (numStates, numStates)))
    for index::Int16 in Int16(1):numStates
        perturbedFreeVariables::Vector{Float64} = copy(X)
        perturbedFreeVariables[index] -= stepSize
        arcMinus::MBD.KArc = propagate(propagator, perturbedFreeVariables, [0, tau], dynamicsModel)
        constraintVectorMinus::StaticArrays.SVector{Int64(numStates), Float64} = StaticArrays.SVector{Int64(numStates), Float64}(getStateByIndex(arcMinus, -1))
        perturbedFreeVariables[index] += 2*stepSize
        arcPlus::MBD.KArc = propagate(propagator, perturbedFreeVariables, [0, tau], dynamicsModel)
        constraintVectorPlus::StaticArrays.SVector{Int64(numStates), Float64} = StaticArrays.SVector{Int64(numStates), Float64}(getStateByIndex(arcPlus, -1))
        STMNumerical[:,index] = (constraintVectorPlus-constraintVectorMinus)./(2*stepSize)
    end
    absDiff::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = STMNumerical.-STMAnalytical
    for r::Int16 in Int16(1):numStates, c::Int16 in Int16(1):numStates
        analytical::Float64 = STMAnalytical[r,c]
        numerical::Float64 = STMNumerical[r,c]
        diff::Float64 = absDiff[r,c]
        useAbs::Bool = ((abs(analytical) < stepSize*1E3) || (abs(numerical) < 1E-12))
        relDiff::Float64 = useAbs ? diff : (diff/abs(numerical))
        errorType::String = useAbs ? "Absolute" : "Relative"
        if relDiff > relTol
            throw(ErrorException("STM error in entry ($r, $c): Expected = $numerical; Actual = $analytical; Difference = $diff; Error = $relDiff ($errorType)"))
        end
    end

    return true
end

"""
    evaluateEquations(dynamicsModel, equationType, t, q)

Return time derivative of state vector

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `equationType::EquationType`: EOM type
- `t::Float64`: Time [ndim]
- `q::Vector{Float64}`: State vector [ndim]
"""
function evaluateEquations(dynamicsModel::KDynamicsModel, equationType::MBD.EquationType, t::Float64, q::Vector{Float64})
    qdot::Vector{Float64} = Vector{Float64}(undef, getStateSize(dynamicsModel, equationType))
    EOMs::MBD.KEquationsOfMotion = getEquationsOfMotion(dynamicsModel, equationType)
    computeDerivatives!(qdot, q, (EOMs,), t)

    return qdot
end

"""
    getCartesianState(dynamicsModel, elementState)

Return Cartesian state from Keplerian orbital elements

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `elementState::Vector{Float64}`: Orbital elements
"""
function getCartesianState(dynamicsModel::KDynamicsModel, elementState::Vector{Float64})
    (a::Float64, e::Float64, i::Float64, Omega::Float64, omega::Float64, theta::Float64) = [elementState...]
    E::Float64 = atan(sqrt(1-e^2)*sin(theta), e+cos(theta))
    r_c::Float64 = a*(1-e*cos(E))
    r_E::Vector{Float64} = [r_c*cos(theta), r_c*sin(theta), 0]
    v_E::Vector{Float64} = (sqrt(getMassParameter(dynamicsModel)*a)/r_c).*[-sin(E), sqrt(1-e^2)*cos(E), 0]
    C::Matrix{Float64} = [cos(omega)*cos(Omega)-sin(omega)*cos(i)*sin(Omega) -sin(omega)*cos(Omega)-cos(omega)*cos(i)*sin(Omega) 0; cos(omega)*sin(Omega)+sin(omega)*cos(i)*cos(Omega) -sin(omega)*sin(Omega)+cos(omega)*cos(i)*cos(Omega) 0; sin(omega)*sin(i) cos(omega)*sin(i) 0]
    
    return [C*r_E; C*v_E]
end

# """
#     getCharLength(dynamicsModel)

# Return CR3BP characteristic length

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# """
# function getCharLength(dynamicsModel::CR3BPDynamicsModel)
#     return getCharLength(dynamicsModel.systemData)
# end

# """
#     getCharTime(dynamicsModel)

# Return CR3BP characteristic time

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# """
# function getCharTime(dynamicsModel::CR3BPDynamicsModel)
#     return getCharTime(dynamicsModel.systemData)
# end

"""
    getEnergy(dynamicsModel, elementState)

Return Keplerian energy

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `elementState::Vector{Float64}`: Orbital elements
"""
function getEnergy(dynamicsModel::KDynamicsModel, elementState::Vector{Float64})
    return -getMassParameter(dynamicsModel)/(2*elementState[1])
end

"""
    getEpochDependencies(dynamicsModel, q)

Return derivative of state with respect to epoch

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getEpochDependencies(dynamicsModel::KDynamicsModel, q_full::Vector{Float64})
    n_full::Int16 = getStateSize(dynamicsModel, MBD.FULL)
    (Int16(length(q_full)) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)

    isEpochIndependent(dynamicsModel) ? (return zeros(Float64, n_simple)) : (return q_full[n_simple*(n_simple+1)+1:n_simple*(n_simple+1)+n_simple])
end

"""
    getEquationsOfMotion(dynamicsModel, equationType)

Return EOMs

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `equationType::EquationType`: EOM type
"""
function getEquationsOfMotion(dynamicsModel::KDynamicsModel, equationType::MBD.EquationType)
    return MBD.KEquationsOfMotion(equationType, dynamicsModel)
end

"""
    getExcursion(dynamicsModel, q)

Return distance from primary

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `q::Vector{Float64}`: State vector [ndim]
"""
function getExcursion(dynamicsModel::KDynamicsModel, q::Vector{Float64})
    primaryPos::Vector{Float64} = getPrimaryState(dynamicsModel)[1:3]

    return LinearAlgebra.norm(q[1:3]-primaryPos)
end

# """
#     getLambertArc(dynamicsModel, initialPos, finalPos, TOF, transferMethod)

# Return initial/final velocities

# # Arguments
# - `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
# - `initalPos::Vector{Float64}`: Initial position [dim]
# - `finalPos::Vector{Float64}`: Final position [dim]
# - `TOF::Float64`: Time of flight [dim]
# - `transferMethod::String`: Transfer method
# """
# function getLambertArc(dynamicsModel::TBPDynamicsModel, initialPos::Vector{Float64}, finalPos::Vector{Float64}, TOF::Float64, transferMethod::String)
#     r0::Float64 = LinearAlgebra.norm(initialPos)
#     rf::Float64 = LinearAlgebra.norm(finalPos)
#     cosdeltanu::Float64 = LinearAlgebra.dot(initialPos, finalPos)/(r0*rf)
#     t_m::Float64 = (transferMethod == "Short") ? 1.0 : -1.0
#     A::Float64 = t_m*sqrt(rf*r0*(1+cosdeltanu))
#     (A == 0) && throw(ErrorException("A = 0 so Lambert arc cannot be computed"))
#     psi_n::Float64 = 0.0
#     c_2::Float64 = 1/2
#     c_3::Float64 = 1/6
#     psi_up::Float64 = 4*pi^2
#     psi_low::Float64 = -8*pi
#     deltat_n::Float64 = TOF+100
#     iter::Int64 = 0
#     while (abs(deltat_n-TOF) >= 1E-6) && (iter < 200)
#         y_n::Float64 = r0+rf+(A*(psi_n*c_3-1))/sqrt(c_2)
#         if (A > 0) && (y_n < 0)
#             psi_low /= 2
#         else
#             x_n::Float64 = sqrt(y_n/c_2)
#             deltat_n = (x_n^3*c_3+A*sqrt(y_n))/sqrt(dynamicsModel.systemData.gravParam)
#             (deltat_n > TOF) ? (psi_up = psi_n) : (psi_low = psi_n)
#         end
#         iter += 1
#         psi_n = (psi_up+psi_low)/2
#         if psi_n > 1E-6
#             c_2 = (1-cos(sqrt(psi_n)))/psi_n
#             c_3 = (sqrt(psi_n)-sin(sqrt(psi_n)))/sqrt(psi_n^3)
#         elseif psi_n < -1E-6
#             c_2 = (1-cosh(sqrt(-psi_n)))/psi_n
#             c_3 = (sinh(sqrt(-psi_n))-sqrt(-psi_n))/sqrt((-psi_n)^3)
#         else
#             c_2 = 1/2
#             c_3 = 1/6
#         end
#     end
#     (iter < 100) || throw(ErrorException("Could not converge Lambert arc"))
#     y_n = r0+rf+(A*(psi_n*c_3-1))/sqrt(c_2)
#     f::Float64 = 1-y_n/r0
#     gdot::Float64 = 1-y_n/rf
#     g::Float64 = A*sqrt(y_n/dynamicsModel.systemData.gravParam)

#     return ((finalPos-f.*initialPos)./g, (gdot.*finalPos-initialPos)./g)
# end

"""
    getMassParameter(dynamicsModel)

Return Keplerian system mass parameter

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
"""
function getMassParameter(dynamicsModel::KDynamicsModel)
    return getMassParameter(dynamicsModel.systemData)
end

"""
    getOrbitalElements(dynamicsModel, state_dim)

Return Keplerian orbital elements from Cartesian state

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `state_dim::Vector{Float64}`: Primary-centered inertial state [dim]
"""
function getOrbitalElements(dynamicsModel::KDynamicsModel, state_dim::Vector{Float64})
    r::Float64 = LinearAlgebra.norm(state_dim[1:3])
    v_r::Float64 = LinearAlgebra.dot(state_dim[4:6], state_dim[1:3])/r
    hvec::Vector{Float64} = LinearAlgebra.cross(state_dim[1:3], state_dim[4:6])
    h::Float64 = LinearAlgebra.norm(hvec)
    i::Float64 = acos(hvec[3]/h)
    n::Vector{Float64} = LinearAlgebra.cross([0, 0, 1], hvec)
    Omega::Float64 = (n[2] < 0) ? 2*pi-acos(n[1]/LinearAlgebra.norm(n)) : acos(n[1]/LinearAlgebra.norm(n))
    evec::Vector{Float64} = LinearAlgebra.cross(state_dim[4:6], hvec)./getMassParameter(dynamicsModel)-state_dim[1:3]./r
    e::Float64 = LinearAlgebra.norm(evec)
    a::Float64 = h^2/(getMassParameter(dynamicsModel)*(1-e^2))
    omega::Float64 = (evec[3] < 0) ? 2*pi-acos(LinearAlgebra.dot(n, evec)/(LinearAlgebra.norm(n)*e)) : acos(LinearAlgebra.dot(n, evec)/(LinearAlgebra.norm(n)*e))
    theta::Float64 = (v_r < 0) ? 2*pi-acos(LinearAlgebra.dot(evec, state_dim[1:3])/(e*r)) : acos(LinearAlgebra.dot(evec, state_dim[1:3])/(e*r))

    return [a, e, i, Omega, omega, theta]
end

"""
    getParameterDependencies(dynamicsModel, q_full)

Return derivative of state with respect to parameters

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getParameterDependencies(dynamicsModel::KDynamicsModel, q_full::Vector{Float64})
    n_full::Int16 = getStateSize(dynamicsModel, MBD.FULL)
    (Int16(length(q_full)) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    i0::Int16 = n_simple*(n_simple+1)
    n::Int16 = n_full-i0
    (n == Int16(0)) && (return zeros(Float64, (n_simple,0)))
    n_params::Int16 = n/n_simple
    dqdp::Matrix{Float64} = zeros(Float64, (n_simple,n_params))
    for r::Int16 in Int16(1):n_simple
        for c::Int16 in Int16(1):n_params
            dqdp[r,c] = q_full[i0+n_simple*(c-1)+r]
        end
    end

    return dqdp
end

"""
    getPeriod(dynamicsModel, elementState)

Return Keplerian orbital period

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `elementState::Vector{Float64}`: Orbital elements
"""
function getPeriod(dynamicsModel::KDynamicsModel, elementState::Vector{Float64})
    return 2*pi/sqrt(getMassParameter(dynamicsModel)/abs(elementState[1])^3)
end

"""
    getPrimaryState(dynamicsModel)

Return state of primary

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
"""
function getPrimaryState(dynamicsModel::KDynamicsModel)
    return zeros(Float64, 6)
end

# """
#     getPsuedopotentialJacobian(dynamicsModel, r)

# Return second derivative of pseudopotential function at given location

# # Arguments
# - `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
# - `r::Vector{Float64}`: Position vector [ndim]
# """
# function getPseudopotentialJacobian(dynamicsModel::KDynamicsModel, r::Vector{Float64})
#     mu::Float64 = getMassRatio(dynamicsModel)
#     r_13::Float64 = sqrt((r[1]+mu)^2+r[2]^2+r[3]^2)
#     r_23::Float64 = sqrt((r[1]-1+mu)^2+r[2]^2+r[3]^2)
#     r_13_3::Float64 = r_13^3
#     r_23_3::Float64 = r_23^3
#     r_13_5::Float64 = r_13_3*r_13^2
#     r_23_5::Float64 = r_23_3*r_23^2
#     ddUdr::Vector{Float64} = zeros(Float64, 6)
#     ddUdr[1] = 1-(1-mu)/r_13_3-mu/r_23_3+3*(1-mu)*(r[1]+mu)^2/r_13_5+3*mu*(r[1]+mu-1)^2/r_23_5
#     ddUdr[2] = 1-(1-mu)/r_13_3-mu/r_23_3+3*(1-mu)*r[2]^2/r_13_5+3*mu*r[2]^2/r_23_5
#     ddUdr[3] = -1*(1-mu)/r_13_3-mu/r_23_3+3*(1-mu)*r[3]^2/r_13_5+3*mu*r[3]^2/r_23_5
#     ddUdr[4] = 3*(1-mu)*(r[1]+mu)*r[2]/r_13_5+3*mu*(r[1]+mu-1)*r[2]/r_23_5
#     ddUdr[5] = 3*(1-mu)*(r[1]+mu)*r[3]/r_13_5+3*mu*(r[1]+mu-1)*r[3]/r_23_5
#     ddUdr[6] = 3*(1-mu)*r[2]*r[3]/r_13_5+3*mu*r[2]*r[3]/r_23_5

#     return ddUdr
# end

"""
    getStateSize(dynamicsModel, equationType)

Return number of state variables

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `equationType::EquationType`: EOM type
"""
function getStateSize(dynamicsModel::KDynamicsModel, equationType::MBD.EquationType)
    type = Dict(MBD.SIMPLE => Int16(6), MBD.STM => Int16(42), MBD.FULL => Int16(42), MBD.ARCLENGTH => Int16(43), MBD.MOMENTUM => Int16(43))

    return type[equationType]
end

"""
    getStateTransitionMatrix(dynamicsModel, q0)

Return STM

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `q0::Vector{Float64}`: Initial state vector with STM in column-major order [ndim]
"""
function getStateTransitionMatrix(dynamicsModel::KDynamicsModel, q0::Vector{Float64})
    n_STM::Int16 = getStateSize(dynamicsModel, MBD.STM)
    (Int16(length(q0)) < n_STM) && throw(ArgumentError("State vector length is $(length(q0)), but should be at least $n_STM"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    STM::Matrix{Float64} = zeros(Float64, (n_simple,n_simple))
    for r::Int16 in Int16(1):n_simple
        for c::Int16 in Int16(1):n_simple
            STM[r,c] = q0[n_simple+n_simple*(c-1)+r]
        end
    end

    return STM
end

"""
    isEpochIndependent(dynamicsModel)

Return true if dynamics model is epoch independent

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
"""
function isEpochIndependent(dynamicsModel::KDynamicsModel)
    return true
end

# """
#     primaryInertialToRotating(dynamicsModel, primary, states_primaryInertial, times)

# Return rotating frame states

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# - `states_primaryInertial::Vector{Vector{Float64}}`: Primary-centered inertial states [ndim]
# - `times::Vector{Float64}`: Epochs [ndim]
# """
# function primaryInertialToRotating(dynamicsModel::CR3BPDynamicsModel, primary::Int64, states_primaryInertial::Vector{Vector{Float64}}, times::Vector{Float64})
#     (length(states_primaryInertial) == length(times)) || throw(ArgumentError("Number of state vectors, $(length(states_primaryInertial)), must match number of times, $(length(times))"))
#     (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
#     states::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
#     for i in 1:length(times)
#         C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([cos(times[i]) -sin(times[i]) 0; sin(times[i]) cos(times[i]) 0; 0 0 1])
#         Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([-sin(times[i]) -cos(times[i]) 0; cos(times[i]) -sin(times[i]) 0; 0 0 0])
#         N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
#         state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(N\states_primaryInertial[i])
#         states[i] = state_primary+getPrimaryState(dynamicsModel, primary)
#     end

#     return states
# end

# """
#     primaryInertial2Rotating(dynamicsModel, secondaryData, states_inertial, times)

# Return rotating frame states

# # Arguments
# - `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
# - `secondaryData::BodyData`: Body data object
# - `states_primaryInertial::Vector{Vector{Float64}}`: Primary-centered inertial states [dim]
# - `times::Vector{Float64}`: Epochs [dim]
# """
# function primaryInertial2Rotating(dynamicsModel::TBPDynamicsModel, secondaryData::MBD.BodyData, states_primaryInertial::Vector{Vector{Float64}}, times::Vector{Float64})
#     (length(states_primaryInertial) == length(times)) || throw(ArgumentError("Number of state vectors, $(length(states_primaryInertial)), must match number of times, $(length(times))"))
#     T::Float64 = 2*pi*secondaryData.orbitRadius^(3/2)/sqrt(dynamicsModel.systemData.gravParam)
#     t::Vector{Float64} = 2*pi*times./T
#     states::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
#     for i in 1:length(times)
#         C::Matrix{Float64} = [cos(t[i]) -sin(t[i]) 0; sin(t[i]) cos(t[i]) 0; 0 0 1]
#         Cdot::Matrix{Float64} = (2*pi/T)*[-sin(t[i]) -cos(t[i]) 0; cos(t[i]) -sin(t[i]) 0; 0 0 0]
#         N::Matrix{Float64} = [C zeros(Float64, (3,3)); Cdot C]
#         states[i] = N\states_primaryInertial[i]
#     end

#     return states
# end

# """
#     rotatingToPrimaryEclipJ2000(dynamicsModel, initialEpoch, states, times)

# Return primary-centered Ecliptic J2000 inertial frame states [ndim]

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `initialEpoch::String`: Initial epoch
# - `states::Vector{Vector{Float64}}`: Rotating states [ndim]
# - `times::Vector{Float64}`: Epochs [ndim]
# """
# function rotatingToPrimaryEclipJ2000(dynamicsModel::CR3BPDynamicsModel, initialEpoch::String, states::Vector{Vector{Float64}}, times::Vector{Float64})
#     numTimes::Int16 = Int16(length(times))
#     (Int16(length(states)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states)), must match number of times, $(length(times))"))
#     lstar::Float64 = getCharLength(dynamicsModel)
#     tstar::Float64 = getCharTime(dynamicsModel)
#     bodyInitialStateDim::Vector{Float64} = getEphemerides(initialEpoch, [0.0], dynamicsModel.systemData.primaryNames[2], dynamicsModel.systemData.primaryNames[1], "ECLIPJ2000")[1][1]
#     primary::MBD.BodyData = dynamicsModel.systemData.primaryData[1]
#     initialEpochTime::Float64 = SPICE.str2et(initialEpoch)
#     bodySPICEElements::StaticArrays.MVector{20, Float64} = StaticArrays.MVector{20, Float64}(SPICE.oscltx(bodyInitialStateDim, initialEpochTime, primary.gravParam))
#     (dynamicsModel.systemData.primaryNames[2] == "Earth") && (bodySPICEElements[3] = 0.0)
#     timesDim::Vector{Float64} = times.*tstar
#     thetadotDim::Float64 = 1/tstar
#     states_primaryInertial::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
#     for i in Int16(1):numTimes
#         state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(states[i]-getPrimaryState(dynamicsModel, 1))
#         state_primaryDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(append!(state_primary[1:3].*lstar, state_primary[4:6].*lstar./tstar))
#         bodyElements::Vector{Float64} = append!([lstar, 0.0], bodySPICEElements[3:5], [bodySPICEElements[6]+timesDim[i]/tstar, initialEpochTime+timesDim[i]], [bodySPICEElements[8]])
#         bodyStateDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(SPICE.conics(bodyElements, initialEpochTime+timesDim[i]))
#         xhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(bodyStateDim[1:3]./lstar)
#         zhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])./LinearAlgebra.norm(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])))
#         yhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(zhat, xhat))
#         C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([xhat yhat zhat])
#         Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([thetadotDim.*yhat -thetadotDim.*xhat zeros(Float64, 3)])
#         N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
#         state_primaryInertialDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(N*state_primaryDim)
#         states_primaryInertial[i] = append!(state_primaryInertialDim[1:3]./lstar, state_primaryInertialDim[4:6].*tstar./lstar)
#     end

#     return states_primaryInertial
# end

# """
#     rotatingToPrimaryInertial(dynamicsModel, primary, states, times)

# Return primary-centered arbitrary inertial frame states [ndim]

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# - `states::Vector{Vector{Float64}}`: Rotating states [ndim]
# - `times::Vector{Float64}`: Epochs [ndim]
# """
# function rotatingToPrimaryInertial(dynamicsModel::CR3BPDynamicsModel, primary::Int64, states::Vector{Vector{Float64}}, times::Vector{Float64})
#     numTimes::Int16 = Int16(length(times))
#     (Int16(length(states)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states)), must match number of times, $(length(times))"))
#     (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
#     states_primaryInertial::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
#     for i in Int16(1):numTimes
#         state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(states[i]-getPrimaryState(dynamicsModel, primary))
#         C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([cos(times[i]-times[1]) -sin(times[i]-times[1]) 0; sin(times[i]-times[1]) cos(times[i]-times[1]) 0; 0 0 1])
#         Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([-sin(times[i]-times[1]) -cos(times[i]-times[1]) 0; cos(times[i]-times[1]) -sin(times[i]-times[1]) 0; 0 0 0])
#         N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
#         states_primaryInertial[i] = N*state_primary
#     end

#     return states_primaryInertial
# end

"""
    solveKeplersEquation(dynamicsModel, elementState)

Return time since periapsis based on true anomaly

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `elementState::Vector{Float64}`: Orbital elements
"""
function solveKeplersEquation(dynamicsModel::KDynamicsModel, elementState::Vector{Float64})
    n::Float64 = sqrt(getMassParameter(dynamicsModel)/abs(elementState[1])^3)
    (elementState[6] == 1.0*pi) && (return pi/n)
    E::Float64 = 2*atan(tan(elementState[6]/2)/sqrt((1+elementState[2])/(1-elementState[2])))
    M::Float64 = E-elementState[2]*sin(E)
    
    M < 0 ? (return (M+2*pi)/n) : (return M/n)
end
