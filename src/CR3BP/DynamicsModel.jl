"""
CR3BPDynamicsModel wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 7/30/23
"""

import MBD: CR3BPDynamicsModel

export appendExtraInitialConditions, evaluateEquations, get2BApproximation
export getEquationsOfMotion, getEpochDependencies, getEquilibriumPoint
export getJacobiConstant, getLinearVariation, getParameterDependencies
export getPrimaryPosition, getPseudopotentialJacobian, getStateSize
export getStateTransitionMatrix, isEpochIndependent, primaryInertial2Rotating

"""
    appendExtraInitialConditions(dynamicsModel, q0_simple, outputEquationType)

Return state vector with extra initial conditions

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `q0_simple::Vector{Float64}`: Simple initial state vector [ndim]
- `outputEquationType::EquationType`: Output state EOM type
"""
function appendExtraInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::MBD.EquationType)
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
    evaluateEquations(dynamicsModel, equationType, t, q, params)

Return time derivative of state vector

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `equationType::EquationType`: EOM type
- `t::Float64`: Time [ndim]
- `q::Vector{Float64}`: State vector [ndim]
- `params::Vector{Float64}`: System parameters
"""
function evaluateEquations(dynamicsModel::CR3BPDynamicsModel, equationType::MBD.EquationType, t::Float64, q::Vector{Float64}, params::Vector{Float64})
    qdot::Vector{Float64} = Vector{Float64}(undef, getStateSize(dynamicsModel, equationType))
    EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, equationType, params)
    computeDerivatives!(qdot, q, EOMs, t)

    return qdot
end

"""
    get2BApproximation(dynamicsModel, bodyData, primary, radius)

Return states of 2BP approximation about primary

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `bodyData::BodyData`: Body data object
- `primary::Int64`: Primary identifier
- `radius::Float64`: Circular radius [ndim]
"""
function get2BApproximation(dynamicsModel::CR3BPDynamicsModel, bodyData::MBD.BodyData, primary::Int64, radius::Float64)
    radius_dim::Float64 = radius*dynamicsModel.systemData.charLength
    circularVelocity_dim::Float64 = sqrt(bodyData.gravParam/radius_dim)
    v::Float64 = circularVelocity_dim*dynamicsModel.systemData.charTime/dynamicsModel.systemData.charLength
    q_primaryInertial::Vector{Float64} = [-radius, 0, 0, 0, v, 0]

    return primaryInertial2Rotating(dynamicsModel, primary, [q_primaryInertial], [0.0])[1]
end

"""
    getEquationsOfMotion(dynamicsModel, equationType; params)

Return EOMs

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `equationType::EquationType`: EOM type
- `params::Vector{Float64}`: System parameters (optional)
"""
function getEquationsOfMotion(dynamicsModel::CR3BPDynamicsModel, equationType::MBD.EquationType, params = nothing)
    return MBD.CR3BPEquationsOfMotion(equationType, dynamicsModel)
end

"""
    getEpochDependencies(dynamicsModel, q)

Return derivative of state with respect to epoch

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getEpochDependencies(dynamicsModel::CR3BPDynamicsModel, q_full::Vector{Float64})
    n_full::Int64 = getStateSize(dynamicsModel, MBD.FULL)
    (length(q_full) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int64 = getStateSize(dynamicsModel, MBD.SIMPLE)

    isEpochIndependent(dynamicsModel) ? (return zeros(Float64, n_simple)) : (return q_full[n_simple*(n_simple+1)+1:n_simple*(n_simple+1)+n_simple])
end

"""
    getEquilibriumPoint(dynamicsModel, point)

Return location of CR3BP equilibirum point in rotating frame

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `point::Int64`: Equilibrium point identifier
"""
function getEquilibriumPoint(dynamicsModel::CR3BPDynamicsModel, point::Int64)
    tol::Float64 = 1E-14
    (1 <= point <= 5) || throw(ArgumentError("Invalid equilibrium point $point"))
    mu::Float64 = getMassRatio(dynamicsModel.systemData)
    pos::Vector{Float64} = zeros(Float64, 3)
    gamma::Float64 = 0.0
    gamma_prev::Float64 = -999
    count::Int64 = 0
    maxCount::Int64 = 20
    if point == 1
        gamma = (mu/(3*(1-mu)))^(1/3)
        while (abs(gamma-gamma_prev) > tol) && (count < maxCount)
            gamma_prev = gamma
            gamma -= (mu/gamma^2-(1-mu)/(1-gamma)^2-gamma-mu+1)/(-2*mu/gamma^3-2*(1-mu)/(1-gamma)^3-1)
            count += 1
        end
        pos[1] = 1-mu-gamma
    elseif point == 2
        gamma = (mu/(3*(1-mu)))^(1/3)
        while (abs(gamma-gamma_prev) > tol) && (count < maxCount)
            gamma_prev = gamma
            gamma -= (-mu/gamma^2-(1-mu)/(1+gamma)^2-mu+1+gamma)/(2*mu/gamma^3+2*(1-mu)/(1+gamma)^3+1)
            count += 1
        end
        pos[1] = 1-mu+gamma
    elseif point == 3
        gamma = 1-7*mu/12
        while (abs(gamma-gamma_prev) > tol) && (count < maxCount)
            gamma_prev = gamma
            gamma -= (mu/(-1-gamma)^2+(1-mu)/gamma^2-mu-gamma)/(-2*mu/(1+gamma)^3-2*(1-mu)/gamma^3-1)
            count += 1
        end
        pos[1] = -mu-gamma
    else
        pos[1] = 1/2-mu
        pos[2] = (point == 4 ? sin(pi/3) : -sin(pi/3))
    end

    (count >= maxCount) ? throw(ErrorException("Could not converge on equilibrium pointlocation")) : (return pos)
end

"""
    getJacobiConstant(dynamicsModel, q)

Return CR3BP Jacobi constant

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `q::Vector{Float64}`: State vector [ndim]
"""
function getJacobiConstant(dynamicsModel::CR3BPDynamicsModel, q::Vector{Float64})
    mu::Float64 = getMassRatio(dynamicsModel.systemData)
    v_2::Float64 = q[4]^2+q[5]^2+q[6]^2
    r_13::Float64 = sqrt((q[1]+mu)^2+q[2]^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-1+mu)^2+q[2]^2+q[3]^2)
    U::Float64 = (1-mu)/r_13+mu/r_23+(1/2)*(q[1]^2+q[2]^2)

    return 2*U-v_2
end

"""
    getLinearVariation(dynamicsModel, equilibriumPos, variation)

Return linear variation about equilibrium point

# Arguments
- `dynamicsModel:::CR3BPDynamicsModel`: CR3BP dynamics model object
- `equilibriumPos::Vector{Float64}`: Equilibrium point position [ndim]
- `variation::Vector{Float64}`: Variation from equilibrium point [ndim]
"""
function getLinearVariation(dynamicsModel::CR3BPDynamicsModel, equilibriumPos::Vector{Float64}, variation::Vector{Float64})
    mu = getMassRatio(dynamicsModel.systemData)
    Uddot::Vector{Float64} = getPseudopotentialJacobian(dynamicsModel, equilibriumPos)
    beta_1::Float64 = 2-(Uddot[1]+Uddot[2])/2
    beta_2_2::Float64 = -Uddot[1]*Uddot[2]
    s::Float64 = sqrt(beta_1+sqrt(beta_1^2+beta_2_2))
    beta_3::Float64 = (s^2+Uddot[1])/(2*s)
    q::Vector{Float64} = push!(equilibriumPos+variation, variation[2]*s/beta_3, -beta_3*variation[1]*s, 0)
    tSpan::Vector{Float64} = [0, 2*pi/s]
    
    return [q, tSpan]
end

"""
    getParameterDependencies(dynamicsModel, q_full)

Return derivative of state with respect to parameters

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getParameterDependencies(dynamicsModel::CR3BPDynamicsModel, q_full::Vector{Float64})
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
    getPrimaryPosition(dynamicsModel, primary)

Return location of primary in rotating frame

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `primary::Int64`: Primary identifier
"""
function getPrimaryPosition(dynamicsModel::CR3BPDynamicsModel, primary::Int64)
    (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
    mu::Float64 = getMassRatio(dynamicsModel.systemData)
    pos::Vector{Float64} = zeros(Float64, 3)
    pos[1] = (primary == 1 ? -mu : (1-mu))

    return pos
end

"""
    getPsuedopotentialJacobian(dynamicsModel, r)

Return second derivative of pseudopotential function at given location

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `r::Vector{Float64}`: Position vector [ndim]
"""
function getPseudopotentialJacobian(dynamicsModel::CR3BPDynamicsModel, r::Vector{Float64})
    mu::Float64 = getMassRatio(dynamicsModel.systemData)
    r_13::Float64 = sqrt((r[1]+mu)^2+r[2]^2+r[3]^2)
    r_23::Float64 = sqrt((r[1]-1+mu)^2+r[2]^2+r[3]^2)
    r_13_3::Float64 = r_13^3
    r_23_3::Float64 = r_23^3
    r_13_5::Float64 = r_13_3*r_13^2
    r_23_5::Float64 = r_23_3*r_23^2
    ddUdr::Vector{Float64} = zeros(Float64, 6)
    ddUdr[1] = 1-(1-mu)/r_13_3-mu/r_23_3+3*(1-mu)*(r[1]+mu)^2/r_13_5+3*mu*(r[1]+mu-1)^2/r_23_5
    ddUdr[2] = 1-(1-mu)/r_13_3-mu/r_23_3+3*(1-mu)*r[2]^2/r_13_5+3*mu*r[2]^2/r_23_5
    ddUdr[3] = -1*(1-mu)/r_13_3-mu/r_23_3+3*(1-mu)*r[3]^2/r_13_5+3*mu*r[3]^2/r_23_5
    ddUdr[4] = 3*(1-mu)*(r[1]+mu)*r[2]/r_13_5+3*mu*(r[1]+mu-1)*r[2]/r_23_5
    ddUdr[5] = 3*(1-mu)*(r[1]+mu)*r[3]/r_13_5+3*mu*(r[1]+mu-1)*r[3]/r_23_5
    ddUdr[6] = 3*(1-mu)*r[2]*r[3]/r_13_5+3*mu*r[2]*r[3]/r_23_5

    return ddUdr
end

"""
    getStateSize(dynamicsModel, equationType)

Return number of state variables

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `equationType::EquationType`: EOM type
"""
function getStateSize(dynamicsModel::CR3BPDynamicsModel, equationType::MBD.EquationType)
    type = Dict(MBD.SIMPLE => 6, MBD.STM => 42, MBD.FULL => 42)

    return type[equationType]
end

"""
    getStateTransitionMatrix(dynamicsModel, q0_STM)

Return STM

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `q0_STM::Vector{Float64}`: Initial state vector with STM in column-major order [ndim]
"""
function getStateTransitionMatrix(dynamicsModel::CR3BPDynamicsModel, q0_STM::Vector{Float64})
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
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
function isEpochIndependent(dynamicsModel::CR3BPDynamicsModel)
    return true
end

"""
    primaryInertial2Rotating(dynamicsModel, primary, states_primaryInertial, times)

Return rotaing frame states

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `primary::Int64`: Primary identifier
- `states_primaryInertial::Vector{Vector{Float64}}`: Primary-centered inertial states [ndim]
- `times::Vector{Float64}`: Epochs [ndim]
"""
function primaryInertial2Rotating(dynamicsModel::CR3BPDynamicsModel, primary::Int64, states_primaryInertial::Vector{Vector{Float64}}, times::Vector{Float64})
    (length(states_primaryInertial) == length(times)) || throw(ArgumentError("Number of state vectors, $(length(states_primaryInertial)), must match number of times, $(length(times))"))
    (1 <= primary <= 2) || throw(ArgumentError("Invalid primary $primary"))
    states::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
    for i in 1:length(times)
        C::Matrix{Float64} = [cos(times[i]) -sin(times[i]) 0; sin(times[i]) cos(times[i]) 0; 0 0 1]
        Cdot::Matrix{Float64} = [-sin(times[i]) -cos(times[i]) 0; cos(times[i]) -sin(times[i]) 0; 0 0 0]
        N::Matrix{Float64} = [C zeros(Float64, (3,3)); Cdot C]
        state_primary::Vector{Float64} = N\states_primaryInertial[i]
        states[i] = state_primary+push!(getPrimaryPosition(dynamicsModel, primary), 0, 0, 0)
    end

    return states
end