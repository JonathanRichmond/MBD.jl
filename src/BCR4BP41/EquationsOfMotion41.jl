"""
BCR4BP P4-B1 equations of motion wrapper

Author: Jonathan Richmond
C: 2/20/25
U: 3/18/25
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
    qdot[7] = sqrt(a4^3/(m4+1))-1
    if params[1].equationType != MBD.SIMPLE
        r_13_5::Float64 = r_13_3*r_13^2
        r_23_5::Float64 = r_23_3*r_23^2
        r_43_5::Float64 = r_43_3*r_43^2
        pseudoPotentialJacobian::StaticArrays.MVector{9, Float64} = StaticArrays.MVector{9, Float64}(zeros(Float64, 9))
        pseudoPotentialJacobian[1] = 1-mu41*omm12/r_13_3-mu14*mu12/r_23_3-omm41/r_43_3+3*mu41*omm12*(q[1]-x1)^2/r_13_5+3*mu41*mu12*(q[1]-x2)^2/r_23_5+3*omm41*(q[1]+mu41)^2/r_43_5
        pseudoPotentialJacobian[2] = 1-mu41*omm12/r_13_3-mu14*mu12/r_23_3-omm41/r_43_3+3*mu41*omm12*(q[2]-y1)^2/r_13_5+3*mu41*mu12*(q[2]-y2)^2/r_23_5+3*omm41*q[2]^2/r_43_5
        pseudoPotentialJacobian[3] = -mu41*omm12/r_13_3-mu14*mu12/r_23_3-omm41/r_43_3+3*mu41*omm12*q[3]^2/r_13_5+3*mu41*mu12*q[3]^2/r_23_5+3*omm41*q[3]^2/r_43_5
        pseudoPotentialJacobian[4] = 3*mu41*omm12*(q[1]-x1)*(q[2]-y1)/r_13_5+3*mu41*mu12*(q[1]-x2)*(q[2]-y2)/r_23_5+3*omm41*(q[1]+mu41)*q[2]/r_43_5
        pseudoPotentialJacobian[5] = 3*mu41*omm12*(q[1]-x1)*q[3]/r_13_5+3*mu41*mu12*(q[1]-x2)*q[3]/r_23_5+3*omm41*(q[1]+mu41)*q[3]/r_43_5
        pseudoPotentialJacobian[6] = 3*mu41*omm12*(q[2]-y1)*q[3]/r_13_5+3*mu41*mu12*(q[2]-y2)*q[3]/r_23_5+3*omm41*q[2]*q[3]/r_43_5
        pseudoPotentialJacobian[7] = mu41*mu12*omm12*sin(q[7])/(a4*r_13_3)-mu41*mu12*omm12*sin(q[7])/(a4*r_23_3)-3*mu41*mu12*omm12*(q[1]-x1)*((q[1]-x1)*sin(q[7])-(q[2]-y1)*cos(q[7]))/(a4*r_13_5)+3*mu14*mu12*omm12*(q[1]-x2)*((q[1]-x2)*sin(q[7])-(q[2]-y2)*cos(q[7]))/(a4*r_23_5)
        pseudoPotentialJacobian[8] = -mu41*mu12*omm12*cos(q[7])/(a4*r_13_3)+mu41*mu12*omm12*cos(q[7])/(a4*r_23_3)-3*mu41*mu12*omm12*(q[2]-y1)*((q[1]-x1)*sin(q[7])-(q[2]-y1)*cos(q[7]))/(a4*r_13_5)+3*mu41*mu12*omm12*(q[2]-y2)*((q[1]-x2)*sin(q[7])-(q[2]-y2)*cos(q[7]))/(a4*r_23_5)
        pseudoPotentialJacobian[9] = -3*mu41*mu12*omm12*q[3]*((q[1]-x1)*sin(q[7])-(q[2]-y1)*cos(q[7]))/(a4*r_13_5)+3*mu41*mu12*omm12*q[3]*((q[1]-x2)*sin(q[7])-(q[2]-y2)*cos(q[7]))/(a4*r_23_5)
        [qdot[7+7*(c-1)+r] = q[10+7*(c-1)+r] for r in 1:3 for c in 1:7]
        for c::Int16 in Int16(1):Int16(7)
            qdot[11+7*(c-1)] = pseudoPotentialJacobian[1]*q[8+7*(c-1)]+pseudoPotentialJacobian[4]*q[9+7*(c-1)]+pseudoPotentialJacobian[5]*q[10+7*(c-1)]+2*q[12+7*(c-1)]+pseudoPotentialJacobian[7]*q[14+7*(c-1)]
            qdot[12+7*(c-1)] = pseudoPotentialJacobian[4]*q[8+7*(c-1)]+pseudoPotentialJacobian[2]*q[9+7*(c-1)]+pseudoPotentialJacobian[6]*q[10+7*(c-1)]-2*q[11+7*(c-1)]+pseudoPotentialJacobian[8]*q[14+7*(c-1)]
            qdot[13+7*(c-1)] = pseudoPotentialJacobian[5]*q[8+7*(c-1)]+pseudoPotentialJacobian[6]*q[9+7*(c-1)]+pseudoPotentialJacobian[3]*q[10+7*(c-1)]+pseudoPotentialJacobian[9]*q[14+7*(c-1)]
        end
        [qdot[7+7*c] = 0 for c in 1:7]
    end
    if params[1].equationType == MBD.ARCLENGTH
        qdot[57] = sqrt(q[4]^2+q[5]^2+q[6]^2)
    elseif params[1].equationType == MBD.MOMENTUM
        qdot[57] = q[1]*q[4]+q[2]*q[5]+q[3]*q[6]
    end
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
