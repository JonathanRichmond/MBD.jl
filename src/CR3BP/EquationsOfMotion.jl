"""
CR3BP equations of motion wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 1/15/25
"""

import StaticArrays
import MBD: CR3BPEquationsOfMotion

export computeDerivatives!, getMassRatio, getStateSize

"""
    computeDerivatives!(qdot, q, params, t)

Return time derivative of state vector

# Arguments
- `qdot::Vector{Float64}`: Time derivative of state vector [ndim]
- `q::Vector{Float64}`: State vector [ndim]
- `params::Tuple{CR3BPEquationsOfMotion, Vararg{Any}}`: Propagation parameters
- `t::Float64`: Time [ndim]
"""
function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, params::Tuple{CR3BPEquationsOfMotion, Vararg{Any}}, t::Float64)
    mu::Float64 = getMassRatio(params[1])
    omm::Float64 = 1-mu
    r_13::Float64 = sqrt((q[1]+mu)^2+q[2]^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-omm)^2+q[2]^2+q[3]^2)
    r_13_3::Float64 = r_13^3
    r_23_3::Float64 = r_23^3
    qdot[1:3] = q[4:6]
    qdot[4] = 2*q[5]+q[1]-omm*(q[1]+mu)/r_13_3-mu*(q[1]-omm)/r_23_3
    qdot[5] = q[2]-2*q[4]-omm*q[2]/r_13_3-mu*q[2]/r_23_3
    qdot[6] = -omm*q[3]/r_13_3-mu*q[3]/r_23_3
    if params[1].equationType != MBD.SIMPLE
        r_13_5::Float64 = r_13_3*r_13^2
        r_23_5::Float64 = r_23_3*r_23^2
        pseudoPotentialJacobian::StaticArrays.MVector{6, Float64} = StaticArrays.MVector{6, Float64}(zeros(Float64, 6))
        pseudoPotentialJacobian[1] = 1-omm/r_13_3-mu/r_23_3+3*omm*(q[1]+mu)^2/r_13_5+3*mu*(q[1]-omm)^2/r_23_5
        pseudoPotentialJacobian[2] = 1-omm/r_13_3-mu/r_23_3+3*omm*q[2]^2/r_13_5+3*mu*q[2]^2/r_23_5
        pseudoPotentialJacobian[3] = -omm/r_13_3-mu/r_23_3+3*omm*q[3]^2/r_13_5+3*mu*q[3]^2/r_23_5
        pseudoPotentialJacobian[4] = 3*omm*(q[1]+mu)*q[2]/r_13_5+3*mu*(q[1]-omm)*q[2]/r_23_5
        pseudoPotentialJacobian[5] = 3*omm*(q[1]+mu)*q[3]/r_13_5+3*mu*(q[1]-omm)*q[3]/r_23_5
        pseudoPotentialJacobian[6] = 3*omm*q[2]*q[3]/r_13_5+3*mu*q[2]*q[3]/r_23_5
        [qdot[6+6*(c-1)+r] = q[9+6*(c-1)+r] for r in 1:3 for c in 1:6]
        for c::Int16 in Int16(1):Int16(6)
            qdot[10+6*(c-1)] = pseudoPotentialJacobian[1]*q[7+6*(c-1)]+pseudoPotentialJacobian[4]*q[8+6*(c-1)]+pseudoPotentialJacobian[5]*q[9+6*(c-1)]+2*q[11+6*(c-1)]
            qdot[11+6*(c-1)] = pseudoPotentialJacobian[4]*q[7+6*(c-1)]+pseudoPotentialJacobian[2]*q[8+6*(c-1)]+pseudoPotentialJacobian[6]*q[9+6*(c-1)]-2*q[10+6*(c-1)]
            qdot[12+6*(c-1)] = pseudoPotentialJacobian[5]*q[7+6*(c-1)]+pseudoPotentialJacobian[6]*q[8+6*(c-1)]+pseudoPotentialJacobian[3]*q[9+6*(c-1)]
        end
    end
    if params[1].equationType == MBD.ARCLENGTH
        qdot[43] = sqrt(q[4]^2+q[5]^2+q[6]^2)
    elseif params[1].equationType == MBD.MOMENTUM
        qdot[43] = q[1]*q[4]+q[2]*q[5]+q[3]*q[6]
    end
end

"""
    getMassRatio(EOMs)

Return CR3BP system mass ratio

# Arguments
- `EOMs::CR3BPEquationsOfMotion`: CR3BP equations of motion object
"""
function getMassRatio(EOMs::CR3BPEquationsOfMotion)
    return getMassRatio(EOMs.dynamicsModel)
end

"""
    getStateSize(EOMs)

Return size of state vector

# Arguments
- `EOMs::CR3BPEquationsOfMotion`: CR3BP equations of motion object
"""
function getStateSize(EOMs::CR3BPEquationsOfMotion)
    return getStateSize(EOMs.dynamicsModel, EOMs.equationType)
end
