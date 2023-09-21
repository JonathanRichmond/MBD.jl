"""
Propagator wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 8/3/23
"""

import DifferentialEquations
import MBD: Propagator

export propagate, propagateWithEvent

"""
    propagate(propagator, q0, tSpan, dynamicsModel; params)

Return propagated arc
# Arguments
- `propagator::Propagator`: Propagator object
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagate(propagator::Propagator, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.AbstractDynamicsModel, params = [])
    arcOut = MBD.Arc(dynamicsModel)
    isempty(params) || setParameters!(arcOut, params)
    EOMs::MBD.AbstractEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType, params)
    q::Vector{Float64} = copy(q0)
    for tIndex::Int64 in 2:length(tSpan)
        if tIndex > 2
            q = copy(getStateByIndex(arcOut, getStateCount(arcOut)))
            deleteStateAndTime!(arcOut, getStateCount(arcOut))
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        arcOut.states = sol.u
        arcOut.times = sol.t
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel; param)

Return propagated arc
# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::ContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `param::Float64`: Propagation parameter
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.ContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.AbstractDynamicsModel, param::Float64)
    arcOut = MBD.Arc(dynamicsModel)
    isempty(params) || setParameters!(arcOut, params)
    EOMs::MBD.AbstractEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType, params)
    q::Vector{Float64} = copy(q0)
    for tIndex::Int64 in 2:length(tSpan)
        if tIndex > 2
            q = copy(getStateByIndex(arcOut, getStateCount(arcOut)))
            deleteStateAndTime!(arcOut, getStateCount(arcOut))
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q, (t0, tf), (EOMs, param))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        arcOut.states = sol.u
        arcOut.times = sol.t
    end

    return arcOut
end
