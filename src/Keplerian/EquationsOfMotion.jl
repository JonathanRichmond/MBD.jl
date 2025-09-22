"""
Keplerian equations of motion wrapper

Author: Jonathan Richmond
C: 9/18/25
"""

import MBD: KEquationsOfMotion

export computeDerivatives!, getMassParameter, getStateSize

"""
    computeDerivatives!(qdot, q, params, t)

Return time derivative of state vector

# Arguments
- `qdot::Vector{Float64}`: Time derivative of state vector [dim]
- `q::Vector{Float64}`: State vector [dim]
- `params::Tuple{KEquationsOfMotion}`: Propagation parameters
- `t::Float64`: Time [dim]
"""
function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, params::Tuple{KEquationsOfMotion}, t::Float64)
    mu::Float64 = getMassParameter(params[1])
    r::Float64 = sqrt(q[1]^2+q[2]^2+q[3]^2)
    qdot[1:3] = q[4:6]
    qdot[4:6] = -mu*q[1:3]./r^3
    if params[1].equationType != MBD.SIMPLE
        A_21::Matrix{Float64} = (3*mu/r^5).*q[1:3]*q[1:3]'-(mu/r^3)*LinearAlgebra.I
        qdot[7:42] = reshape([zeros(Float64, (3,3)) LinearAlgebra.I; A_21 zeros(Float64, (3,3))], 36)
    end
    if params[1].equationType == MBD.ARCLENGTH
        qdot[43] = sqrt(q[4]^2+q[5]^2+q[6]^2)
    elseif params[1].equationType == MBD.MOMENTUM
        qdot[43] = q[1]*q[4]+q[2]*q[5]+q[3]*q[6]
    end
end

"""
    getMassParameter(EOMs)

Return Keplerian system mass parameter

# Arguments
- `EOMs::KEquationsOfMotion`: Keplerian equations of motion object
"""
function getMassParameter(EOMs::KEquationsOfMotion)
    return getMassParameter(EOMs.dynamicsModel)
end

"""
    getStateSize(EOMs)

Return size of state vector

# Arguments
- `EOMs::KEquationsOfMotion`: Keplerian equations of motion object
"""
function getStateSize(EOMs::KEquationsOfMotion)
    return getStateSize(EOMs.dynamicsModel, EOMs.equationType)
end
