"""
Event functions

Author: Jonathan Richmond
C: 9/20/23
U: 7/16/24
"""

import DifferentialEquations

export arclengthCondition, momentumDifferenceCondition, p1DistanceCondition
export p2DistanceCondition, terminateAffect!, xzPlaneCrossingCondition
export zValueCondition

"""
    arclengthCondition(state, time, integrator)

Return event condition for specified arclength

# Arguments
- `state::Vector{Float64}`: State vector with arclength [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [arclength]
"""
function arclengthCondition(state::Vector{Float64}, time::Float64, integrator)
    state[43]-integrator.p[2]
end

"""
    momentumDifferenceCondition(state, time, integrator)

Return event condition for specified momentum integral difference

# Arguments
- `state::Vector{Float64}`: State vector with momentum integral [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [propagator, dynamicsModel, q0, momentumDifference]
"""
function momentumDifferenceCondition(state::Vector{Float64}, time::Float64, integrator)
    orbitArc::MBD.Arc = propagate(integrator.p[2], appendExtraInitialConditions(integrator.p[3], integrator.p[4], MBD.MOMENTUM), [0, time], integrator.p[3])
    abs(state[43]-getStateByIndex(orbitArc, -1)[43])-integrator.p[5]
end

"""
    p1DistanceCondition(state, time, integrator)

Return event condition for specified distance from P1

# Arguments
- `state::Vector{Float64}`: State vector [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [p1Distance]
"""
function p1DistanceCondition(state::Vector{Float64}, time::Float64, integrator)
    mu::Float64 = integrator.p[1].mu
    sqrt((state[1]+mu)^2+state[2]^2+state[3]^2)-integrator.p[2]
end

"""
    p2DistanceCondition(state, time, integrator)

Return event condition for specified distance from P2

# Arguments
- `state::Vector{Float64}`: State vector [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [p2Distance]
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
- `integrator`: Integrator object with params: [z]
"""
function zValueCondition(state::Vector{Float64}, time::Float64, integrator)
    state[3]-integrator.p[2]
end
