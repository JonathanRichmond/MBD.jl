"""
TBP equations of motion wrapper

Author: Jonathan Richmond
C: 9/14/23
U: 5/15/24
"""

import LinearAlgebra
import MBD: TBPEquationsOfMotion

export computeDerivatives!

"""
    computeDerivatives!(qdot, q, params, t)

Return time derivative of state vector

# Arguments
- `qdot::Vector{Float64}`: Time derivative of state vector [ndim]
- `q::Vector{Float64}`: State vector [ndim]
- `params::Tuple{TBPEquationsOfMotion}`: Propagation parameters
- `t::Float64`: Time [ndim]
"""
function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, params::Tuple{TBPEquationsOfMotion}, t::Float64)
    r::Float64 = sqrt(q[1]^2+q[2]^2+q[3]^2)
    qdot[1:3] = q[4:6]
    qdot[4:6] = -1*q[1:3]./r^3
    if (params[1].equationType == MBD.STM) || (params[1].equationType == MBD.FULL) || (params[1].equationType == MBD.ARCLENGTH)
        A_21::Matrix{Float64} = (3/r^5).*q[1:3]*q[1:3]'-(1/r^3)*LinearAlgebra.I
        qdot[7:42] = reshape([zeros(Float64, (3,3)) LinearAlgebra.I; A_21 zeros(Float64, (3,3))]', (1,36))
    end
    if params[1].equationType == MBD.ARCLENGTH
        qdot[43] = sqrt(q[4]^2+q[5]^2+q[6]^2)
    end
end
