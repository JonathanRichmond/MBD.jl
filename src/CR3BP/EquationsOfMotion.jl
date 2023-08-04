"""
CR3BPEquationsOfMotion wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 7/30/23
"""

import MBD: CR3BPEquationsOfMotion

export computeDerivatives!

"""
    computeDerivatives!(EOM, qdot, q, t)

Return time derivative of state vector

# Arguments
- `EOMs::CR3BPEquationsOfMotion`: EOM obbject
- `qdot::Vector{Float64}`: Time derivative of state vector [ndim]
- `q::Vector{Float64}`: State vector [ndim]
- `t::Float64`: Time [ndim]
"""
function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, EOMs::CR3BPEquationsOfMotion, t::Float64)
    omm::Float64 = 1-EOMs.mu
    r_13::Float64 = sqrt((q[1]+EOMs.mu)^2+q[2]^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-omm)^2+q[2]^2+q[3]^2)
    r_13_3::Float64 = r_13^3
    r_23_3::Float64 = r_23^3
    qdot[1:3] = q[4:6]
    qdot[4] = 2*q[5]+q[1]-omm*(q[1]+EOMs.mu)/r_13_3-EOMs.mu*(q[1]-omm)/r_23_3
    qdot[5] = q[2]-2*q[1]-omm*q[2]/r_13_3-EOMs.mu*q[2]/r_23_3
    qdot[6] = -omm*q[3]/r_13_3-EOMs.mu*q[3]/r_23_3
    if (EOMs.equationType == MBD.STM) || (EOMs.equationType == MBD.FULL)
        r_13_5::Float64 = r_13_3*r_13^2
        r_23_5::Float64 = r_23_3*r_23^2
        pseudoPotentialJacobian::Vector{Float64} = Vector{Float64}(undef, 6)
        pseudoPotentialJacobian[1] = 1-omm/r_13_3-EOMs.mu/r_23_3+3*omm*(q[1]+EOMs.mu)^2/r_13_5+3*EOMs.mu*(q[1]-omm)^2/r_23_5
        pseudoPotentialJacobian[2] = 1-omm/r_13_3-EOMs.mu/r_23_3+3*omm*q[2]^2/r_13_5+3*EOMs.mu*q[2]^2/r_23_5
        pseudoPotentialJacobian[3] = -omm/r_13_3-EOMs.mu/r_23_3+3*omm*q[3]^2/r_13_5+3*EOMs.mu*q[3]^2/r_23_5
        pseudoPotentialJacobian[4] = 3*omm*(q[1]+EOMs.mu)*q[2]/r_13_5+3*EOMs.mu*(q[1]-omm)*q[2]/r_23_5
        pseudoPotentialJacobian[5] = 3*omm*(q[1]+EOMs.mu)*q[3]/r_13_5+3*EOMs.mu*(q[1]-omm)*q[3]/r_23_5
        pseudoPotentialJacobian[6] = 3*omm*q[2]*q[3]/r_13_5+3*EOMs.mu*q[2]*q[3]/r_23_5
        [qdot[6+6*(c-1)+r] = q[9+6*(c-1)+r] for r in 1:3 for c in 1:6]
        for c in 1:6
            qdot[10+6*(c-1)] = pseudoPotentialJacobian[1]*q[7+6*(c-1)]+pseudoPotentialJacobian[4]*q[8+6*(c-1)]+pseudoPotentialJacobian[5]*q[9+6*(c-1)]+2*q[11+6*(c-1)]
            qdot[11+6*(c-1)] = pseudoPotentialJacobian[4]*q[7+6*(c-1)]+pseudoPotentialJacobian[2]*q[8+6*(c-1)]+pseudoPotentialJacobian[6]*q[9+6*(c-1)]-2*q[10+6*(c-1)]
            qdot[12+6*(c-1)] = pseudoPotentialJacobian[5]*q[7+6*(c-1)]+pseudoPotentialJacobian[6]*q[8+6*(c-1)]+pseudoPotentialJacobian[3]*q[9+6*(c-1)]
        end
    end
end
