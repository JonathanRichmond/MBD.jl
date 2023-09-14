"""
TBP equations of motion wrapper

Author: Jonathan Richmond
C: 9/14/23
"""

import LinearAlgebra
import MBD: TBPEquationsOfMotion

export computeDerivatives!

"""
    computeDerivatives!(qdot, q, EOMs, t)

Return time derivative of state vector

# Arguments
- `qdot::Vector{Float64}`: Time derivative of state vector [ndim]
- `q::Vector{Float64}`: State vector [ndim]
- `EOMs::TBPEquationsOfMotion`: TBP EOM obbject
- `t::Float64`: Time [ndim]
"""
function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, EOMs::TBPEquationsOfMotion, t::Float64)
    r::Float64 = sqrt(q[1]^2+q[2]^2+q[3]^2)
    qdot[1:3] = q[4:6]
    qdot[4:6] = -1*q[1:3]./r^3
    if (EOMs.equationType == MBD.STM) || (EOMs.equationType == MBD.FULL)
        A_21 = (3/r^5).*q[1:3]*q[1:3]'-(1/r^3)*LinearAlgebra.I
        qdot[7:42] = vec([zeros(Float64, (3,3)) LinearAlgebra.I; A_21 zeros(Float64, (3,3))])
    end
end
