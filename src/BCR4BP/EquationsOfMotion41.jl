"""
BCR4BP P4-B1 equations of motion wrapper

Author: Jonathan Richmond
C: 2/20/25
"""

import MBD: BCR4BP41EquationsOfMotion

export computeDerivatives!, getStateSize, get12MassRatio, get4Distance, get4Mass, get41MassRatio

"""
    computeDerivatives!(qdot, q, params, t)

Return time derivative of state vector

# Arguments
- `qdot::Vector{Float64}`: Time derivative of state vector [ndim]
- `q::Vector{Float64}`: State vector [ndim]
- `params::Tuple{BCR4BP41EquationsOfMotion, Vararg{Any}}`: Propagation parameters
- `t::Float64`: Time [ndim]
"""
function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, params::Tuple{BCR4BP41EquationsOfMotion, Vararg{Any}}, t::Float64)
    mu12::Float64 = get12MassRatio(params[1])
    mu41::Float64 = get41MassRatio(params[1])
    omm12::Float64 = 1-mu12
    omm41::Float64 = 1-mu41
    a4::Float64 = get4Distance(params[1])
    m4::Float64 = get4Mass(params[1])
    x1::Float64 = omm41-mu12*cos(q[7])/a4
    y1::Float64 = -mu12*sin(q[7])/a4
    x2::Float64 = omm41+omm12*cos(q[7])/a4
    y2::Float64 = omm12*sin(q[7])/a4
    r_13::Float64 = sqrt((q[1]-x1)^2+(q[2]-y1)^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-x2)^2+(q[2]-y2)^2+q[3]^2)
    r_43::Float64 = sqrt((q[1]+mu41)^2+q[2]^2+q[3]^2)
    r_13_3::Float64 = r_13^3
    r_23_3::Float64 = r_23^3
    r_43_3::Float64 = r_43^3
    qdot[1:3] = q[4:6]
    qdot[4] = 2*q[5]+q[1]-mu41*omm12*(q[1]-x1)/r_13_3-mu41*mu12*(q[1]-x2)/r_23_3-omm41*(q[1]+mu41)/r_43_3
    qdot[5] = q[2]-2*q[4]-mu41*omm12*(q[2]-y1)/r_13_3-mu41*mu12*(q[2]-y2)/r_23_3-omm41*q[2]/r_43_3
    qdot[6] = -mu41*omm12*q[3]/r_13_3-mu41*mu12*q[3]/r_23_3-omm41*q[3]/r_43_3
    qdot[7] = 1-sqrt((m4+1)/(a4^3))
#     if params[1].equationType != MBD.SIMPLE
#         r_13_5::Float64 = r_13_3*r_13^2
#         r_23_5::Float64 = r_23_3*r_23^2
#         pseudoPotentialJacobian::StaticArrays.MVector{6, Float64} = StaticArrays.MVector{6, Float64}(zeros(Float64, 6))
#         pseudoPotentialJacobian[1] = 1-omm/r_13_3-mu/r_23_3+3*omm*(q[1]+mu)^2/r_13_5+3*mu*(q[1]-omm)^2/r_23_5
#         pseudoPotentialJacobian[2] = 1-omm/r_13_3-mu/r_23_3+3*omm*q[2]^2/r_13_5+3*mu*q[2]^2/r_23_5
#         pseudoPotentialJacobian[3] = -omm/r_13_3-mu/r_23_3+3*omm*q[3]^2/r_13_5+3*mu*q[3]^2/r_23_5
#         pseudoPotentialJacobian[4] = 3*omm*(q[1]+mu)*q[2]/r_13_5+3*mu*(q[1]-omm)*q[2]/r_23_5
#         pseudoPotentialJacobian[5] = 3*omm*(q[1]+mu)*q[3]/r_13_5+3*mu*(q[1]-omm)*q[3]/r_23_5
#         pseudoPotentialJacobian[6] = 3*omm*q[2]*q[3]/r_13_5+3*mu*q[2]*q[3]/r_23_5
#         [qdot[6+6*(c-1)+r] = q[9+6*(c-1)+r] for r in 1:3 for c in 1:6]
#         for c::Int16 in Int16(1):Int16(6)
#             qdot[10+6*(c-1)] = pseudoPotentialJacobian[1]*q[7+6*(c-1)]+pseudoPotentialJacobian[4]*q[8+6*(c-1)]+pseudoPotentialJacobian[5]*q[9+6*(c-1)]+2*q[11+6*(c-1)]
#             qdot[11+6*(c-1)] = pseudoPotentialJacobian[4]*q[7+6*(c-1)]+pseudoPotentialJacobian[2]*q[8+6*(c-1)]+pseudoPotentialJacobian[6]*q[9+6*(c-1)]-2*q[10+6*(c-1)]
#             qdot[12+6*(c-1)] = pseudoPotentialJacobian[5]*q[7+6*(c-1)]+pseudoPotentialJacobian[6]*q[8+6*(c-1)]+pseudoPotentialJacobian[3]*q[9+6*(c-1)]
#         end
#     end
#     if params[1].equationType == MBD.ARCLENGTH
#         qdot[43] = sqrt(q[4]^2+q[5]^2+q[6]^2)
#     elseif params[1].equationType == MBD.MOMENTUM
#         qdot[43] = q[1]*q[4]+q[2]*q[5]+q[3]*q[6]
#     end
end

"""
    getStateSize(EOMs)

Return size of state vector

# Arguments
- `EOMs::BCR4BP41EquationsOfMotion`: BCR4BP P4-B1 equations of motion object
"""
function getStateSize(EOMs::BCR4BP41EquationsOfMotion)
    return getStateSize(EOMs.dynamicsModel, EOMs.equationType)
end

"""
    get12MassRatio(EOMs)

Return BCR4BP P1-P2 system mass ratio

# Arguments
- `EOMs::BCR4BP41EquationsOfMotion`: BCR4BP P4-B1 equations of motion object
"""
function get12MassRatio(EOMs::BCR4BP41EquationsOfMotion)
    return get12MassRatio(EOMs.dynamicsModel)
end

"""
    get4Distance(EOMs)

Return BCR4BP P4 distance from B1 [ndim]

# Arguments
- `EOMs::BCR4BP41EquationsOfMotion`: BCR4BP P4-B1 equations of motion object
"""
function get4Distance(EOMs::BCR4BP41EquationsOfMotion)
    return get4Distance(EOMs.dynamicsModel)
end

"""
    get4Mass(EOMs)

Return BCR4BP P4 mass [ndim]

# Arguments
- `EOMs::BCR4BP41EquationsOfMotion`: BCR4BP P4-B1 equations of motion object
"""
function get4Mass(EOMs::BCR4BP41EquationsOfMotion)
    return get4Mass(EOMs.dynamicsModel)
end

"""
    get41MassRatio(EOMs)

Return BCR4BP P4-B1 system mass ratio

# Arguments
- `EOMs::BCR4BP41EquationsOfMotion`: BCR4BP P4-B1 equations of motion object
"""
function get41MassRatio(EOMs::BCR4BP41EquationsOfMotion)
    return get41MassRatio(EOMs.dynamicsModel)
end
