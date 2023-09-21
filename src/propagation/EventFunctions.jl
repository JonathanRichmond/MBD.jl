"""
Event functions

Author: Jonathan Richmond
C: 9/20/23
"""

import DifferentialEquations

export arclengthCondition!, terminateAffect!

"""
    arclengthCondition!(state, time, integrator)

Return event condition

# Arguments
- `state::Vector{Float64}`: State vector with arclength [ndim]
- `time::Float64`: Time [ndim]
- `integrator::DEIntegrator`: Integrator object
"""
function arclengthCondition!(state::Vector{Float64}, time::Float64, integrator::DifferentialEquations.DEIntegrator)
    state[43]-integrator.opts.userdata[:conditionValue]
end

"""
    terminateAffect!(integrator)

Return event effect

# Arguments
- `integrator::DEIntegrator`: Integrator object
"""
function terminateAffect!(integrator::DifferentialEquations.DEIntegrator)
    DifferentialEquations.terminate!(integrator)
end
