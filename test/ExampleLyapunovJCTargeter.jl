"""
Jacobi constant perpendicular crossing targeter for CR3BP Lyapunov orbits

Author: Jonathan Richmond
C: 1/11/23
U: 9/10/23
"""

using MBD, LinearAlgebra

export LyapunovJCTargeter, correct, propagateState

"""
    LyapunovJCTargeter(dynamicsModel)

CR3BP Lyapunov Jacobi constant targeter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
mutable struct LyapunovJCTargeter <: MBD.AbstractTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                   # CR3BP dynamics model

    function LyapunovJCTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        return new(dynamicsModel)
    end
end

"""
    correct(targeter, q0, tSpan, target)

Return corrected multiple shooter problem

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Lyapunov Jacobi constant targeter object
- `q0::Vector{Float64}`: Initial conditions [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `targetJC::Float64`: Target Jacobi constant
"""
function correct(targeter::LyapunovJCTargeter, q0::Vector{Float64}, tSpan::Vector{Float64}, targetJC::Float64)
    halfPeriod::Float64 = (tSpan[2]-tSpan[1])/2
    qf::Vector{Float64} = propagateState(targeter, q0, [0.0, halfPeriod])
    originNode = MBD.Node(tSpan[1], q0, targeter.dynamicsModel)
    originNode.state.name = "Initial State"
    terminalNode = MBD.Node(halfPeriod, qf, targeter.dynamicsModel)
    terminalNode.state.name = "Target State"
    segment = MBD.Segment(halfPeriod, originNode, terminalNode)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    problem = MBD.MultipleShooterProblem()
    addSegment!(problem, segment)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    JCConstraint = MBD.JacobiConstraint(originNode, targetJC)
    qfConstraint = MBD.StateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, qfConstraint)
    addConstraint!(problem, continuityConstraint)
    shooter = MBD.MultipleShooter()
    solved::MBD.MultipleShooterProblem = solve!(shooter, problem)
    convergenceCheck::MBD.ConstraintVectorL2NormConvergenceCheck = shooter.convergenceCheck
    for constraint::MBD.AbstractConstraint in getConstraints(solved)
        (LinearAlgebra.norm(evaluateConstraint(constraint, getFreeVariableIndexMap!(solved), getFreeVariableVector!(solved))) > convergenceCheck.maxVectorNorm) && println("ERROR: Multiple shooter failed to converge for $(typeof(constraint))")
    end

    return solved
end

"""
    propagateState(targeter, q0, tSpan)

Return propagated state

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Lyapunov Jacobi constant targeter object
- `q0::Vector{Float64}`: Initial conditions [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::LyapunovJCTargeter, q0::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.Arc = propagate(propagator, q0, tSpan, targeter.dynamicsModel)

    return copy(getStateByIndex(arc, getStateCount(arc)))
end
