"""
Event functions

Author: Jonathan Richmond
C: 9/20/23
U: 2/15/24
"""

import DifferentialEquations

export arclengthCondition, p2DistanceCondition, terminateAffect!
export xzPlaneCrossingCondition, zValueCondition

"""
    arclengthCondition(state, time, integrator)

Return event condition for specified arclength

# Arguments
- `state::Vector{Float64}`: State vector with arclength [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object
"""
function arclengthCondition(state::Vector{Float64}, time::Float64, integrator)
    state[43]-integrator.p[2]
end

"""
    p2DistanceCondition(state, time, integrator)

Return event condition for specified distance

# Arguments
- `state::Vector{Float64}`: State vector [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object
"""
function p2DistanceCondition(state::Vector{Float64}, time::Float64, integrator)
    mu::Float64 = integrator.p[1].mu
    sqrt((state[1]-1+mu)^2+state[2]^2+state[3]^2)-integrator.p[2]
end

"""
    terminateAffect!(integrator)

Return event effect of termination

# Arguments
- `integrator`: Integrator object
"""
function terminateAffect!(integrator)
    DifferentialEquations.terminate!(integrator)
end

"""
    xzPlaneCrossingCondition(state, time, integrator)

Return event condition for xz-plane crossing

# Arguments
- `state::Vector{Float64}`: State vector [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object
"""
function xzPlaneCrossingCondition(state::Vector{Float64}, time::Float64, integrator)
    state[2]
end

"""
    zValueCondition(state, time, integrator)

Return event condition for z-value crossing

# Arguments
- `state::Vector{Float64}`: State vector [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object
"""
function zValueCondition(state::Vector{Float64}, time::Float64, integrator)
    state[3]-integrator.p[2]
end
