"""
Propagator wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 10/27/25
"""

import DifferentialEquations
import MBD: Propagator

export propagate, propagateWithEvent, propagateWithEvents, propagateWithPeriodicEvent

"""
    propagate(propagator, q0, tSpan, dynamicsModel)

Return propagated Keplerian arc

# Arguments
- `propagator::Propagator`: Propagator object
- `q0::Vector{Float64}`: Initial state vector [dim]
- `tSpan::Vector{Float64}`: Time span [dim]
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
"""
function propagate(propagator::Propagator, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.KDynamicsModel)
    arcOut = MBD.KArc(dynamicsModel)
    EOMs::MBD.KEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

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
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
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
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagate(propagator, q0, tSpan, dynamicsModel)

Return propagated BCR4BP P4-B1 arc

# Arguments
- `propagator::Propagator`: Propagator object
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function propagate(propagator::Propagator, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP41DynamicsModel)
    arcOut = MBD.BCR4BP41Arc(dynamicsModel)
    EOMs::MBD.BCR4BP41EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs,))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::ContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.ContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.CR3BPDynamicsModel, params = [])
    arcOut = MBD.CR3BPArc(dynamicsModel)
    EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::ContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.ContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP12DynamicsModel, params = [])
    arcOut = MBD.BCR4BP12Arc(dynamicsModel)
    EOMs::MBD.BCR4BP12EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::ContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.ContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP41DynamicsModel, params = [])
    arcOut = MBD.BCR4BP41Arc(dynamicsModel)
    EOMs::MBD.BCR4BP41EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::VectorContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.VectorContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.CR3BPDynamicsModel, params = [])
    arcOut = MBD.CR3BPArc(dynamicsModel)
    EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithEvent(propagator, callbackEvent, q0, tSpan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::VectorContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithEvent(propagator::Propagator, callbackEvent::DifferentialEquations.VectorContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP12DynamicsModel, params = [])
    arcOut = MBD.BCR4BP12Arc(dynamicsModel)
    EOMs::MBD.BCR4BP12EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithEvents(propagator, callbackEvent, q0, tSpan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::VectorContinuousCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `params::Vector{Any}`: Propagation parameters (optional)
"""
function propagateWithEvents(propagator::Propagator, callbackEvent::DifferentialEquations.VectorContinuousCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP41DynamicsModel, params = [])
    arcOut = MBD.BCR4BP41Arc(dynamicsModel)
    EOMs::MBD.BCR4BP41EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, copy(q0), (t0, tf), (EOMs, :none, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return (arcOut, sol.prob.p[2])
end

"""
    propagateWithPeriodicEvent(propagator, callbackEvent, q0, tspan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::DiscreteCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithPeriodicEvent(propagator::Propagator, callbackEvent::DifferentialEquations.DiscreteCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.CR3BPDynamicsModel, params = [])
    arcOut = MBD.CR3BPArc(dynamicsModel)
    EOMs::MBD.CR3BPEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithPeriodicEvent(propagator, callbackEvent, q0, tspan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::DiscreteCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithPeriodicEvent(propagator::Propagator, callbackEvent::DifferentialEquations.DiscreteCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP12DynamicsModel, params = [])
    arcOut = MBD.BCR4BP12Arc(dynamicsModel)
    EOMs::MBD.BCR4BP12EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end

"""
    propagateWithPeriodicEvent(propagator, callbackEvent, q0, tspan, dynamicsModel, params)

Return propagated arc

# Arguments
- `propagator::Propagator`: Propagator object
- `callbackEvent::DiscreteCallback`: Propagation callback
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `params::Vector{Float64}`: Propagation parameters (optional)
"""
function propagateWithPeriodicEvent(propagator::Propagator, callbackEvent::DifferentialEquations.DiscreteCallback, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.BCR4BP41DynamicsModel, params = [])
    arcOut = MBD.BCR4BP41Arc(dynamicsModel)
    EOMs::MBD.BCR4BP41EquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType)
    for tIndex::Int16 in Int16(2):Int16(length(tSpan))
        if tIndex > Int16(2)
            q0 = copy(getStateByIndex(arcOut, -1))
            deleteStateAndTime!(arcOut, -1)
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem = DifferentialEquations.ODEProblem(computeDerivatives!, q0, (t0, tf), (EOMs, params...))
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, callback = callbackEvent, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        append!(arcOut.states, sol.u)
        append!(arcOut.times, sol.t)
    end

    return arcOut
end
