"""
Propagator wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 1/15/25
"""

import DifferentialEquations
import MBD: Propagator

export propagate, propagateWithEvent

"""
    propagate(propagator, q0, tSpan, dynamicsModel)

Return propagated CR3BP arc

# Arguments
- `propagator::Propagator`: Propagator object
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
function propagate(propagator::Propagator, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.CR3BPDynamicsModel)
    arcOut = MBD.CR3BPArc(dynamicsModel)
    EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, getStateCount(arcOut)))
            deleteStateAndTime!(arcOut, getStateCount(arcOut))
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        arcOut.states = sol.u
        arcOut.times = sol.t
    end

    return arcOut
end

"""
    propagate(propagator, q0, tSpan, dynamicsModel)

Return propagated BCR4BP P1-P2 arc

# Arguments
- `propagator::Propagator`: Propagator object
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
function propagate(propagator::Propagator, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP12DynamicsModel)
    arcOut = MBD.BCR4BP12Arc(dynamicsModel)
    EOMs::MBD.BCR4BP12EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, getStateCount(arcOut)))
            deleteStateAndTime!(arcOut, getStateCount(arcOut))
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        arcOut.states = sol.u
        arcOut.times = sol.t
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::ContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.ContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.CR3BPDynamicsModel, params = [])
    arcOut = MBD.CR3BPArc(dynamicsModel)
    EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, getStateCount(arcOut)))
            deleteStateAndTime!(arcOut, getStateCount(arcOut))
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        arcOut.states = sol.u
        arcOut.times = sol.t
    end

    return arcOut
end
