"""
Propagator wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 8/3/23
"""

import DifferentialEquations
import MBD: Propagator

export propagate

"""
    propagate(propagator, q0, tSpan, dynamicsModel; params)

Return propagated arc
# Arguments
- `propagator::Propagator`: Propagator object
- `q0::Vector{Float64}`: Initial state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `params::Vector{Float64}`: System parameters (optional)
"""
function propagate(propagator::Propagator, q0::Vector{Float64}, tSpan::Vector{Float64}, dynamicsModel::MBD.AbstractDynamicsModel, params = [])
    arcOut = MBD.Arc(dynamicsModel)
    #isempty(params) || setParameters!(arcOut, params)
    EOMs::MBD.AbstractEquationsOfMotion = getEquationsOfMotion(dynamicsModel, propagator.equationType, params)
    q::Vector{Float64} = copy(q0)
    for tIndex::Int64 in 2:length(tSpan)
        if tIndex > 2
            #q = copy(getStateByIndex(arcOut, getStateCount(arcOut)))
            #deleteStateAndTime!(arcOut, getStateCount(arcOut))
        end
        t0::Float64 = tSpan[tIndex-1]
        tf::Float64 = tSpan[tIndex]
        problem::DifferentialEquations.ODEProblem = DifferentialEquations.ODEProblem(computeDerivatives!, q, (t0, tf), EOMs)
        sol::DifferentialEquations.ODESolution = DifferentialEquations.solve(problem, propagator.integratorFactory.integrator, abstol = propagator.absTol, reltol = propagator.relTol, dtmax = propagator.maxStep, maxiters = propagator.maxEvaluationCount)
        arcOut.states = sol.u
        arcOut.times = sol.t
    end

    return arcOut
end
