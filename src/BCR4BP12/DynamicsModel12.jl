"""
BCR4BP P1-P2 dynamics model wrapper

Author: Jonathan Richmond
C: 2/26/25
U: 3/3/25
"""

import StaticArrays
import MBD: BCR4BP12DynamicsModel

export appendExtraInitialConditions, checkSTM, getEquationsOfMotion, getStateSize, get12MassRatio
export get4Distance, get4Mass, rotating122Rotating41

"""
    appendExtraInitialConditions(dynamicsModel, q0_simple, outputEquationType)

Return state vector with extra initial conditions

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `q0_simple::Vector{Float64}`: Simple initial state vector [ndim]
- `outputEquationType::EquationType`: Output state EOM type
"""
function appendExtraInitialConditions(dynamicsModel::BCR4BP12DynamicsModel, q0_simple::Vector{Float64}, outputEquationType::MBD.EquationType)
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
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `relTol::Float64`: Relative tolerance (default = 1E-5)
"""
function checkSTM(dynamicsModel::BCR4BP12DynamicsModel, relTol::Float64 = 1E-5)
    stepSize::Float64 = sqrt(eps(Float64))
    numStates::Int16 = Int16(getStateSize(dynamicsModel, MBD.SIMPLE))
    propagator = MBD.Propagator()
    propagatorSTM = MBD.Propagator(equationType = MBD.STM)
    X::Vector{Float64} = [0.9, 0, 0, 0, -0.7, 0, 0]
    tau::Float64 = 0.1
    arc::MBD.BCR4BP12Arc = propagate(propagatorSTM, appendExtraInitialConditions(dynamicsModel, X, MBD.STM), [0, tau], dynamicsModel)
    STMAnalytical::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64}(reshape(getStateByIndex(arc, -1)[(numStates+1):(numStates+numStates^2)], (7,7)))
    STMNumerical::StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64}(zeros(Float64, (numStates, numStates)))
    for index::Int16 in Int16(1):numStates
        perturbedFreeVariables::Vector{Float64} = copy(X)
        perturbedFreeVariables[index] -= stepSize
        arcMinus::MBD.BCR4BP12Arc = propagate(propagator, perturbedFreeVariables, [0, tau], dynamicsModel)
        constraintVectorMinus::StaticArrays.SVector{Int64(numStates), Float64} = StaticArrays.SVector{Int64(numStates), Float64}(getStateByIndex(arcMinus, -1))
        perturbedFreeVariables[index] += 2*stepSize
        arcPlus::MBD.BCR4BP12Arc = propagate(propagator, perturbedFreeVariables, [0, tau], dynamicsModel)
        constraintVectorPlus::StaticArrays.SVector{Int64(numStates), Float64} = StaticArrays.SVector{Int64(numStates), Float64}(getStateByIndex(arcPlus, -1))
        STMNumerical[:,index] = (constraintVectorPlus-constraintVectorMinus)./(2*stepSize)
    end
    absDiff::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64}(STMNumerical-STMAnalytical)
    relDiff::StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64}(copy(absDiff))
    for r::Int16 in Int16(1):numStates
        for c::Int16 in Int16(1):numStates
            if abs(STMAnalytical[r,c]) < stepSize*1E3
                relDiff[r,c] = absDiff[r,c]
            elseif abs(STMNumerical[r,c]) > 1E-11
                relDiff[r,c] = absDiff[r,c]/STMNumerical[r,c]
            end
            if abs(relDiff[r,c]) > relTol
                throw(ErrorException("Jacobian error in entry ($r, $c): Expected = $(jacobianNumerical[r,c]); Actual = $(jacobianAnalytical[r,c]) (Relative error = $(relDiff[r,c]))"))
                return false
            end
        end
    end

    return true
end

# """
#     evaluateEquations(dynamicsModel, equationType, t, q)

# Return time derivative of state vector

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `equationType::EquationType`: EOM type
# - `t::Float64`: Time [ndim]
# - `q::Vector{Float64}`: State vector [ndim]
# """
# function evaluateEquations(dynamicsModel::CR3BPDynamicsModel, equationType::MBD.EquationType, t::Float64, q::Vector{Float64})
#     qdot::Vector{Float64} = Vector{Float64}(undef, getStateSize(dynamicsModel, equationType))
#     EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, equationType)
#     computeDerivatives!(qdot, q, (EOMs,), t)

#     return qdot
# end

# """
#     get2BApproximation(dynamicsModel, bodyData, primary, radius)

# Return states of 2BP approximation about primary

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `bodyData::BodyData`: Body data object
# - `primary::Int64`: Primary identifier
# - `radius::Float64`: Circular radius [ndim]
# """
# function get2BApproximation(dynamicsModel::CR3BPDynamicsModel, bodyData::MBD.BodyData, primary::Int64, radius::Float64)
#     lstar::Float64 = getCharLength(dynamicsModel)
#     tstar::Float64 = getCharTime(dynamicsModel)
#     radius_dim::Float64 = radius*lstar
#     circularVelocity_dim::Float64 = sqrt(bodyData.gravParam/radius_dim)
#     v::Float64 = circularVelocity_dim*tstar/lstar
#     q_primaryInertial::Vector{Float64} = [-radius, 0, 0, 0, v, 0]

#     return primaryInertial2Rotating(dynamicsModel, primary, [q_primaryInertial], [0.0])[1]
# end

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

# """
#     getEpochDependencies(dynamicsModel, q)

# Return derivative of state with respect to epoch

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `q_full::Vector{Float64}`: Full state vector [ndim]
# """
# function getEpochDependencies(dynamicsModel::CR3BPDynamicsModel, q_full::Vector{Float64})
#     n_full::Int16 = getStateSize(dynamicsModel, MBD.FULL)
#     (Int64(length(q_full)) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
#     n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)

#     isEpochIndependent(dynamicsModel) ? (return zeros(Float64, n_simple)) : (return q_full[n_simple*(n_simple+1)+1:n_simple*(n_simple+1)+n_simple])
# end

"""
    getEquationsOfMotion(dynamicsModel, equationType)

Return EOMs

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `equationType::EquationType`: EOM type
"""
function getEquationsOfMotion(dynamicsModel::BCR4BP12DynamicsModel, equationType::MBD.EquationType)
    return MBD.BCR4BP12EquationsOfMotion(equationType, dynamicsModel)
end

# """
#     getEquilibriumPoint(dynamicsModel, point)

# Return location of CR3BP equilibirum point in rotating frame

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `point::Int64`: Equilibrium point identifier
# """
# function getEquilibriumPoint(dynamicsModel::CR3BPDynamicsModel, point::Int64)
#     tol::Float64 = 1E-14
#     (1 <= point <= 5) || throw(ArgumentError("Invalid equilibrium point $point"))
#     mu::Float64 = getMassRatio(dynamicsModel)
#     pos::StaticArrays.MVector{3, Float64} = StaticArrays.MVector{3}(zeros(Float64, 3))
#     gamma::Float64 = 0.0
#     gamma_prev::Float64 = -999.0
#     count::Int16 = 0
#     maxCount::Int16 = 20
#     if point == 1
#         gamma = (mu/(3*(1-mu)))^(1/3)
#         while (abs(gamma-gamma_prev) > tol) && (count < maxCount)
#             gamma_prev = gamma
#             gamma -= (mu/gamma^2-(1-mu)/(1-gamma)^2-gamma-mu+1)/(-2*mu/gamma^3-2*(1-mu)/(1-gamma)^3-1)
#             count += 1
#         end
#         pos[1] = 1-mu-gamma
#     elseif point == 2
#         gamma = (mu/(3*(1-mu)))^(1/3)
#         while (abs(gamma-gamma_prev) > tol) && (count < maxCount)
#             gamma_prev = gamma
#             gamma -= (-mu/gamma^2-(1-mu)/(1+gamma)^2-mu+1+gamma)/(2*mu/gamma^3+2*(1-mu)/(1+gamma)^3+1)
#             count += 1
#         end
#         pos[1] = 1-mu+gamma
#     elseif point == 3
#         gamma = 1-7*mu/12
#         while (abs(gamma-gamma_prev) > tol) && (count < maxCount)
#             gamma_prev = gamma
#             gamma -= (mu/(-1-gamma)^2+(1-mu)/gamma^2-mu-gamma)/(-2*mu/(1+gamma)^3-2*(1-mu)/gamma^3-1)
#             count += 1
#         end
#         pos[1] = -mu-gamma
#     else
#         pos[1] = 1/2-mu
#         pos[2] = (point == 4 ? sin(pi/3) : -sin(pi/3))
#     end

#     (count >= maxCount) ? throw(ErrorException("Could not converge on equilibrium pointlocation")) : (return pos)
# end

# """
#     getExcursion(dynamicsModel, primary, q)

# Return distance from primary

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# - `q::Vector{Float64}`: State vector [ndim]
# """
# function getExcursion(dynamicsModel::CR3BPDynamicsModel, primary::Int64, q::Vector{Float64})
#     lstar::Float64 = getCharLength(dynamicsModel)
#     primaryPos::Vector{Float64} = getPrimaryPosition(dynamicsModel, primary)

#     return LinearAlgebra.norm(q[1:3]-primaryPos)*lstar
# end

# """
#     getJacobiConstant(dynamicsModel, q)

# Return CR3BP Jacobi constant

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `q::Vector{Float64}`: State vector [ndim]
# """
# function getJacobiConstant(dynamicsModel::CR3BPDynamicsModel, q::Vector{Float64})
#     mu::Float64 = getMassRatio(dynamicsModel)
#     v_2::Float64 = q[4]^2+q[5]^2+q[6]^2
#     r_13::Float64 = sqrt((q[1]+mu)^2+q[2]^2+q[3]^2)
#     r_23::Float64 = sqrt((q[1]-1+mu)^2+q[2]^2+q[3]^2)
#     U::Float64 = (1-mu)/r_13+mu/r_23+(1/2)*(q[1]^2+q[2]^2)

#     return 2*U-v_2
# end

# """
#     getParameterDependencies(dynamicsModel, q_full)

# Return derivative of state with respect to parameters

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `q_full::Vector{Float64}`: Full state vector [ndim]
# """
# function getParameterDependencies(dynamicsModel::CR3BPDynamicsModel, q_full::Vector{Float64})
#     n_full::Int16 = getStateSize(dynamicsModel, MBD.FULL)
#     (Int16(length(q_full)) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
#     n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
#     i0::Int16 = n_simple*(n_simple+1+(isEpochIndependent(dynamicsModel) ? 0 : 1))
#     n::Int16 = n_full-i0
#     (n == Int16(0)) && (return zeros(Float64, (n_simple,0)))
#     n_params::Int16 = n/n_simple
#     dqdp::Matrix{Float64} = zeros(Float64, (n_simple,n_params))
#     for r::Int16 in Int16(1):n_simple
#         for c::Int16 in Int16(1):n_params
#             dqdp[r,c] = q_full[i0+n_simple*(c-1)+r]
#         end
#     end

#     return dqdp
# end

# """
#     getPrimaryPosition(dynamicsModel, primary)

# Return location of primary in rotating frame

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# """
# function getPrimaryPosition(dynamicsModel::CR3BPDynamicsModel, primary::Int64)
#     (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
#     mu::Float64 = getMassRatio(dynamicsModel)
#     pos::Vector{Float64} = zeros(Float64, 3)
#     pos[1] = (primary == 1 ? -mu : (1-mu))

#     return pos
# end

# """
#     getPsuedopotentialJacobian(dynamicsModel, r)

# Return second derivative of pseudopotential function at given location

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `r::Vector{Float64}`: Position vector [ndim]
# """
# function getPseudopotentialJacobian(dynamicsModel::CR3BPDynamicsModel, r::Vector{Float64})
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
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `equationType::EquationType`: EOM type
"""
function getStateSize(dynamicsModel::BCR4BP12DynamicsModel, equationType::MBD.EquationType)
    type = Dict(MBD.SIMPLE => Int16(7), MBD.STM => Int16(56), MBD.FULL => Int16(56), MBD.ARCLENGTH => Int16(57), MBD.MOMENTUM => Int16(57))

    return type[equationType]
end

# """
#     getStateTransitionMatrix(dynamicsModel, q0)

# Return STM

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `q0::Vector{Float64}`: Initial state vector with STM in column-major order [ndim]
# """
# function getStateTransitionMatrix(dynamicsModel::CR3BPDynamicsModel, q0::Vector{Float64})
#     n_STM::Int16 = getStateSize(dynamicsModel, MBD.STM)
#     (Int16(length(q0)) < n_STM) && throw(ArgumentError("State vector length is $(length(q0)), but should be at least $n_STM"))
#     n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
#     STM::Matrix{Float64} = zeros(Float64, (n_simple,n_simple))
#     for r::Int16 in Int16(1):n_simple
#         for c::Int16 in Int16(1):n_simple
#             STM[r,c] = q0[n_simple+n_simple*(c-1)+r]
#         end
#     end

#     return STM
# end

# """
#     getTidalAcceleration(dynamicsModel, primary, r)

# Return tidal acceleration due to primary at given location

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# - `r::Vector{Float64}`: Position vector [ndim]
# """
# function getTidalAcceleration(dynamicsModel::CR3BPDynamicsModel, primary::Int64, r::Vector{Float64})
#     (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
#     mu::Float64 = getMassRatio(dynamicsModel)
#     r_13::Float64 = sqrt((r[1]+mu)^2+r[2]^2+r[3]^2)
#     r_23::Float64 = sqrt((r[1]-1+mu)^2+r[2]^2+r[3]^2)
#     aMultiplier::Float64 = (primary == 1) ? -(1-mu)/r_13^3 : -mu/r_23^3
#     xDist::Float64 = (primary == 1) ? r[1]+mu : r[1]-1+mu

#     return [aMultiplier*xDist+r[1], aMultiplier*r[2]+r[2], aMultiplier*r[3]]
# end

"""
    get12MassRatio(dynamicsModel)

Return BCR4BP P1-P2 system mass ratio

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
function get12MassRatio(dynamicsModel::BCR4BP12DynamicsModel)
    return get12MassRatio(dynamicsModel.systemData)
end

"""
    get4Distance(dynamicsModel)

Return BCR4BP P4 distance from B1 [ndim]

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
function get4Distance(dynamicsModel::BCR4BP12DynamicsModel)
    return get4Distance(dynamicsModel.systemData)
end

"""
    get4Mass(dynamicsModel)

Return BCR4BP P4 mass [ndim]

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
function get4Mass(dynamicsModel::BCR4BP12DynamicsModel)
    return get4Mass(dynamicsModel.systemData)
end

# """
#     isEpochIndependent(dynamicsModel)

# Return true if dynamics model is epoch independent

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# """
# function isEpochIndependent(dynamicsModel::CR3BPDynamicsModel)
#     return true
# end

# """
#     primaryInertial2Rotating(dynamicsModel, primary, states_primaryInertial, times)

# Return rotating frame states

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# - `states_primaryInertial::Vector{Vector{Float64}}`: Primary-centered inertial states [ndim]
# - `times::Vector{Float64}`: Epochs [ndim]
# """
# function primaryInertial2Rotating(dynamicsModel::CR3BPDynamicsModel, primary::Int64, states_primaryInertial::Vector{Vector{Float64}}, times::Vector{Float64})
#     (length(states_primaryInertial) == length(times)) || throw(ArgumentError("Number of state vectors, $(length(states_primaryInertial)), must match number of times, $(length(times))"))
#     (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
#     states::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
#     for i in 1:length(times)
#         C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3}([cos(times[i]) -sin(times[i]) 0; sin(times[i]) cos(times[i]) 0; 0 0 1])
#         Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3}([-sin(times[i]) -cos(times[i]) 0; cos(times[i]) -sin(times[i]) 0; 0 0 0])
#         N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6}([C zeros(Float64, (3,3)); Cdot C])
#         state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6}(N\states_primaryInertial[i])
#         states[i] = state_primary+push!(getPrimaryPosition(dynamicsModel, primary), 0, 0, 0)
#     end

#     return states
# end

"""
    rotating122Rotating41(dynamicsModel, states12, times12)

Return BCR4BP P4-B1 rotating frame states and times [ndim]

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `states12::Vector{Vector{Float64}}`: BCR4BP P1-P2 rotating frame states [ndim]
- `times12::Vector{Float64}`: BCR4BP P1-P2 rotating frame times [ndim]
"""
function rotating122Rotating41(dynamicsModel::BCR4BP12DynamicsModel, states12::Vector{Vector{Float64}}, times12::Vector{Float64})
    numTimes::Int16 = Int16(length(times12))
    m4::Float64 = get4Mass(dynamicsModel.systemData)
    a4::Float64 = get4Distance(dynamicsModel.systemData)
    theta4dot::Float64 = sqrt((m4+1)/(a4^3))-1
    states41::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
    for t::Int16 = Int16(1):numTimes
        state::StaticArrays.SVector{7, Float64} = StaticArrays.SVector{7, Float64}(states12[t])
        theta4::Float64 = state[7]
        C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([-cos(theta4) -sin(theta4) 0; sin(theta4) -cos(theta4) 0; 0 0 1])
        Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}(theta4dot.*[sin(theta4) -cos(theta4) 0; cos(theta4) sin(theta4) 0; 0 0 0])
        states41[t] = [(1/a4).*C zeros(Float64, 3, 4); sqrt(a4/(m4+1)).*Cdot sqrt(a4/(m4+1)).*C zeros(Float64, 3, 1); zeros(Float64, 1, 6) -1]*state+append!([1-get41MassRatio(dynamicsModel.systemData)], zeros(Float64, 5), [pi])
    end
    times41::Vector{Float64} = (get12CharTime(dynamicsModel.systemData)/get41CharTime(dynamicsModel.systemData)).*times12

    return (states41, times41)
end

# """
#     rotating2PrimaryEclipJ2000(dynamicsModel, initialEpoch, states, times)

# Return primary-centered Ecliptic J2000 inertial frame states [ndim]

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `initialEpoch::String`: Initial epoch
# - `states::Vector{Vector{Float64}}`: Rotating states [ndim]
# - `times::Vector{Float64}`: Epochs [ndim]
# """
# function rotating2PrimaryEclipJ2000(dynamicsModel::CR3BPDynamicsModel, initialEpoch::String, states::Vector{Vector{Float64}}, times::Vector{Float64})
#     numTimes::Int16 = Int16(length(times))
#     (Int16(length(states)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states)), must match number of times, $(length(times))"))
#     lstar::Float64 = getCharLength(dynamicsModel)
#     tstar::Float64 = getCharTime(dynamicsModel)
#     bodyInitialStateDim::Vector{Vector{Float64}} = getEphemerides(initialEpoch, [0.0], dynamicsModel.systemData.primaryNames[2], dynamicsModel.systemData.primaryNames[1], "ECLIPJ2000")[1]
#     primary = MBD.BodyData(dynamicsModel.systemData.primaryNames[1])
#     initialEpochTime::Float64 = SPICE.str2et(initialEpoch)
#     bodySPICEElements::StaticArrays.SVector{20, Float64} = StaticArrays.SVector{20, Float64}(SPICE.oscltx(bodyInitialStateDim[1], initialEpochTime, primary.gravParam))
#     (dynamicsModel.systemData.primaryNames[2] == "Earth") && (bodySPICEElements[3] = 0.0)
#     timesDim::Vector{Float64} = times.*tstar
#     thetadotDim::Float64 = 1/tstar
#     states_primaryInertial::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
#     for i in Int16(1):numTimes
#         stateDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(append!(states[i][1:3].*lstar, states[i][4:6].*lstar./tstar))
#         state_primaryDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(stateDim-push!(getPrimaryPosition(dynamicsModel, 1).*lstar, 0, 0, 0))
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
#     rotating2PrimaryInertial(dynamicsModel, primary, states, times)

# Return primary-centered arbitrary inertial frame states [ndim]

# # Arguments
# - `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
# - `primary::Int64`: Primary identifier
# - `states::Vector{Vector{Float64}}`: Rotating states [ndim]
# - `times::Vector{Float64}`: Epochs [ndim]
# """
# function rotating2PrimaryInertial(dynamicsModel::CR3BPDynamicsModel, primary::Int64, states::Vector{Vector{Float64}}, times::Vector{Float64})
#     numTimes::Int16 = Int16(length(times))
#     (Int16(length(states)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states)), must match number of times, $(length(times))"))
#     (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
#     states_primaryInertial::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
#     for i in Int16(1):numTimes
#         state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(states[i]-push!(getPrimaryPosition(dynamicsModel, primary), 0, 0, 0))
#         C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([cos(times[i]) -sin(times[i]) 0; sin(times[i]) cos(times[i]) 0; 0 0 1])
#         Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([-sin(times[i]) -cos(times[i]) 0; cos(times[i]) -sin(times[i]) 0; 0 0 0])
#         N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
#         states_primaryInertial[i] = N*state_primary
#     end

#     return states_primaryInertial
# end
