"""
Event functions

Author: Jonathan Richmond
C: 9/20/23
"""

import DifferentialEquations

export arclengthCondition, terminateAffect!

"""
    arclengthCondition(state, time, integrator)

Return event condition

# Arguments
- `state::Vector{Float64}`: State vector with arclength [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object
"""
function arclengthCondition(state::Vector{Float64}, time::Float64, integrator::DifferentialEquations.Integrator)
    state[43]-integrator.p[2]
end

"""
    terminateAffect!(integrator)

Return event effect

# Arguments
- `integrator`: Integrator object
"""
function terminateAffect!(integrator)
    DifferentialEquations.terminate!(integrator)
end
