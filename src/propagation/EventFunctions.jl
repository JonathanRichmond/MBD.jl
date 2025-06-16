"""
Event functions

Author: Jonathan Richmond
C: 9/20/23
U: 6/16/25
"""

import DifferentialEquations, LinearAlgebra, Logging

export arclengthCondition, momentumDifferenceConditionBCR4BP12, momentumDifferenceConditionBCR4BP41
export momentumDifferenceConditionCR3BP, primaryDistanceCondition, p1DistanceCondition
export p2DistanceCondition, renormalize!, terminateAffect!, xzPlaneCrossingCondition
export zValueCondition

"""
    arclengthCondition(state, time, integrator)

Return event condition for specified arclength

# Arguments
- `state::Vector{Float64}`: State vector with arclength [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [dynamicsModel, arclength]
"""
function arclengthCondition(state::Vector{Float64}, time::Float64, integrator)
    n_arclength::Int16 = getStateSize(integrator.p[2], MBD.ARCLENGTH)
    state[n_arclength]-integrator.p[3]
end

"""
    momentumDifferenceConditionBCR4BP12(state, time, integrator)

Return event condition for specified momentum integral difference

# Arguments
- `state::Vector{Float64}`: State vector with momentum integral [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [propagator, dynamicsModel, q0, momentumDifference]
"""
function momentumDifferenceConditionBCR4BP12(state::Vector{Float64}, time::Float64, integrator)
    n_momentum::Int16 = getStateSize(integrator.p[3], MBD.MOMENTUM)
    orbitArc::MBD.BCR4BP12Arc = propagate(integrator.p[2], appendExtraInitialConditions(integrator.p[3], integrator.p[4], MBD.MOMENTUM), [0, time], integrator.p[3])
    abs(state[n_momentum]-getStateByIndex(orbitArc, -1)[n_momentum])-integrator.p[5]
end

"""
    momentumDifferenceConditionBCR4BP41(state, time, integrator)

Return event condition for specified momentum integral difference

# Arguments
- `state::Vector{Float64}`: State vector with momentum integral [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [propagator, dynamicsModel, q0, momentumDifference]
"""
function momentumDifferenceConditionBCR4BP41(state::Vector{Float64}, time::Float64, integrator)
    n_momentum::Int16 = getStateSize(integrator.p[3], MBD.MOMENTUM)
    orbitArc::MBD.BCR4BP41Arc = propagate(integrator.p[2], appendExtraInitialConditions(integrator.p[3], integrator.p[4], MBD.MOMENTUM), [0, time], integrator.p[3])
    abs(state[n_momentum]-getStateByIndex(orbitArc, -1)[n_momentum])-integrator.p[5]
end

"""
    momentumDifferenceConditionCR3BP(state, time, integrator)

Return event condition for specified momentum integral difference

# Arguments
- `state::Vector{Float64}`: State vector with momentum integral [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [propagator, dynamicsModel, q0, momentumDifference]
"""
function momentumDifferenceConditionCR3BP(state::Vector{Float64}, time::Float64, integrator)
    n_momentum::Int16 = getStateSize(integrator.p[3], MBD.MOMENTUM)
    orbitArc::MBD.CR3BPArc = propagate(integrator.p[2], appendExtraInitialConditions(integrator.p[3], integrator.p[4], MBD.MOMENTUM), [0, time], integrator.p[3])
    abs(state[n_momentum]-getStateByIndex(orbitArc, -1)[n_momentum])-integrator.p[5]
end

"""
    primaryDistanceCondition(output, state, time, integrator)

Return event conditions for specified distances from primaries

# Arguments
- `output`: Condition output vector []
- `state::Vector{Float64}`: State vector [ndim]
- `time::Float64`: Time [ndim]
- `integrator`: Integrator object with params: [dynamicsModel, p1Distance, p2Distance, p4Distance]
"""
function primaryDistanceCondition(output, state::Vector{Float64}, time::Float64, integrator)
    r1::Vector{Float64} = getPrimaryState(integrator.p[2], 1, state[7])[1:3]
    r2::Vector{Float64} = getPrimaryState(integrator.p[2], 2, state[7])[1:3]
    r4::Vector{Float64} = getPrimaryState(integrator.p[2], 4, state[7])[1:3]
    d1::Float64 = LinearAlgebra.norm(state[1:3]-r1)
    d2::Float64 = LinearAlgebra.norm(state[1:3]-r2)
    d4::Float64 = LinearAlgebra.norm(state[1:3]-r4)
    output[1:3] = [d1-integrator.p[3], d2-integrator.p[4], d4-integrator.p[5]]
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
    mu::Float64 = getMassRatio(integrator.p[1])
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
    mu::Float64 = getMassRatio(integrator.p[1])
    sqrt((state[1]-1+mu)^2+state[2]^2+state[3]^2)-integrator.p[2]
end

"""
    renormalize!(integrator)

Return event effect of STM renormalization

# Arguments
- `integrator`: Integrator object with params: [dynamicsModel, Rs]
"""
function renormalize!(integrator)
    n_simple::Int16 = getStateSize(integrator.p[2], MBD.SIMPLE)
    n_STM::Int16 = getStateSize(integrator.p[2], MBD.STM)
    Phi::Matrix{Float64} = reshape(integrator.u[(n_simple+1):n_STM], (n_simple,n_simple))
    F = LinearAlgebra.qr(Phi)
    integrator.u[(n_simple+1):n_STM] = vec(Matrix(F.Q))
    push!(integrator.p[3], Matrix(F.R))
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
    terminateAffect!(integrator, index)

Return event effect of termination

# Arguments
- `integrator`: Integrator object with params: [dynamicsModel, ...]
- `index::Int64`: Condition index
"""
function terminateAffectIndex!(integrator, index)
    Logging.@info "Propagation terminated with crash into $(integrator.p[2].systemData.primaryNames[index])"
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
