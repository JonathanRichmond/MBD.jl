"""
Initial x-state perpendicular crossing targeter for CR3BP Lyapunov orbits

Author: Jonathan Richmond
C: 1/27/25
"""

using MBD, LinearAlgebra, StaticArrays

export LyapunovJCTargeter
export correct, getMonodromy, getPeriod, getPeriodicOrbit, propagateState

"""
    Lyapunovx0Targeter(dynamicsModel)

CR3BP Lyapunov initial x-state targeter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct Lyapunovx0Targeter <: MBD.AbstractTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model

    function Lyapunovx0Targeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, q0, tSpan, targetx0; tol)

Return corrected multiple shooter problem

# Arguments
- `targeter::Lyapunovx0Targeter`: CR3BP Lyapunov initial x-state targeter object
- `q0::Vector{Float64}`: Initial conditions [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `targetx0::Float64`: Target initial x-state [ndim]
- `tol::Float64`: Convergence tolerance (optional)
"""
function correct(targeter::Lyapunovx0Targeter, q0::Vector{Float64}, tSpan::Vector{Float64}, targetx0::Float64, tol::Float64 = 1E-11)
    halfPeriod::Float64 = (tSpan[2]-tSpan[1])/2
    qf::Vector{Float64} = propagateState(targeter, q0, [0, halfPeriod])
    originNode = MBD.CR3BPNode(tSpan[1], q0, targeter.dynamicsModel)
    originNode.state.name = "Initial State"
    terminalNode = MBD.CR3BPNode(halfPeriod, qf, targeter.dynamicsModel)
    terminalNode.state.name = "Target State"
    segment = MBD.CR3BPSegment(halfPeriod, originNode, terminalNode)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    problem = MBD.CR3BPMultipleShooterProblem()
    addSegment!(problem, segment)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    x0Constraint = MBD.CR3BPStateConstraint(originNode, [1], [targetx0])
    qfConstraint = MBD.CR3BPStateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(problem, x0Constraint)
    addConstraint!(problem, qfConstraint)
    addConstraint!(problem, continuityConstraint)
    shooter = MBD.CR3BPMultipleShooter(tol)
    # shooter.printProgress = true
    solved::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)
    convergenceCheck::MBD.ConstraintVectorL2NormConvergenceCheck = shooter.convergenceCheck
    for constraint::MBD.AbstractConstraint in getConstraints(solved)
        (LinearAlgebra.norm(evaluateConstraint(constraint, getFreeVariableIndexMap!(solved), getFreeVariableVector!(solved))) > convergenceCheck.maxVectorNorm) && println("ERROR: Multiple shooter failed to converge for $(typeof(constraint))")
    end

    return solved
end

"""
    getMonodromy(targeter, problem)

Return orbit monodromy matrix

# Arguments
- `targeter::Lyapunovx0Targeter`: CR3BP Lyapunov initial x-state targeter object
- `problem::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getMonodromy(targeter::Lyapunovx0Targeter, problem::MBD.CR3BPMultipleShooterProblem)
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    initialGuess::Vector{Float64} = appendExtraInitialConditions(targeter.dynamicsModel, problem.nodes[1].state.data, MBD.STM)
    arc::MBD.CR3BPArc = propagate(propagator, initialGuess, [0, problem.segments[1].TOF.data[1]], targeter.dynamicsModel)
    halfPeriodState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(arc, -1))
    halfPeriodSTM::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([halfPeriodState[7:12] halfPeriodState[13:18] halfPeriodState[19:24] halfPeriodState[25:30] halfPeriodState[31:36] halfPeriodState[37:42]])
    G::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([1.0 0 0 0 0 0; 0 -1.0 0 0 0 0; 0 0 1.0 0 0 0; 0 0 0 -1.0 0 0; 0 0 0 0 1.0 0; 0 0 0 0 0 -1.0])
    Omega::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([0 1.0 0; -1.0 0 0; 0 0 0])
    A::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([zeros(Float64, (3,3)) [-1.0 0 0; 0 -1.0 0; 0 0 -1.0]; [1.0 0 0; 0 1.0 0; 0 0 1.0] -2*Omega])
    B::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([-2*Omega [1.0 0 0; 0 1.0 0; 0 0 1.0]; [-1.0 0 0; 0 -1.0 0; 0 0 -1.0] zeros(Float64, (3,3))])
    
    return G*A*(halfPeriodSTM')*B*G*halfPeriodSTM
end

"""
    getPeriod(targeter, problem)

Return orbit period

# Arguments
- `targeter::Lyapunovx0Targeter`: CR3BP Lyapunov initial x-state targeter object
- `problem::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getPeriod(targeter::Lyapunovx0Targeter, problem::MBD.CR3BPMultipleShooterProblem)
    return 2*problem.segments[1].TOF.data[1]
end

"""
    getPeriodicOrbit(targeter, family, orbit)

Return periodic orbit object

# Arguments
- `targeter::Lyapunovx0Targeter`: CR3BP Lyapunov initial x-state targeter object
- `family::CR3BPContinuationFamily`: CR3BP continuation family object
- `orbit::Int64`: Orbit identifier
"""
function getPeriodicOrbit(targeter::Lyapunovx0Targeter, family::MBD.CR3BPContinuationFamily, orbit::Int64)
    period::Float64 = 2*family.segments[orbit][1].TOF.data[1]
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    initialGuess::Vector{Float64} = appendExtraInitialConditions(targeter.dynamicsModel, family.nodes[orbit][1].state.data, MBD.STM)
    arc::MBD.CR3BPArc = propagate(propagator, initialGuess, [0, family.segments[orbit][1].TOF.data[1]], targeter.dynamicsModel)
    halfPeriodState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(arc, -1))
    halfPeriodSTM::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([halfPeriodState[7:12] halfPeriodState[13:18] halfPeriodState[19:24] halfPeriodState[25:30] halfPeriodState[31:36] halfPeriodState[37:42]])
    G::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([1.0 0 0 0 0 0; 0 -1.0 0 0 0 0; 0 0 1.0 0 0 0; 0 0 0 -1.0 0 0; 0 0 0 0 1.0 0; 0 0 0 0 0 -1.0])
    Omega::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([0 1.0 0; -1.0 0 0; 0 0 0])
    A::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([zeros(Float64, (3,3)) [-1.0 0 0; 0 -1.0 0; 0 0 -1.0]; [1.0 0 0; 0 1.0 0; 0 0 1.0] -2*Omega])
    B::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([-2*Omega [1.0 0 0; 0 1.0 0; 0 0 1.0]; [-1.0 0 0; 0 -1.0 0; 0 0 -1.0] zeros(Float64, (3,3))])
    monodromy::Matrix{Float64} = G*A*(halfPeriodSTM')*B*G*halfPeriodSTM

    return MBD.CR3BPPeriodicOrbit(targeter.dynamicsModel, initialGuess[1:6], period, monodromy)
end

# """
#     getTangentBifurcationOrbit(targeter, orbitFamily, l, h)

# Return tangent bifurcation orbit

# # Arguments
# - `targeter::LyapunovJCTargeter`: CR3BP Laypunov Jacobi constant targeter object
# - `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
# - `l::Int64`: Low orbit identifier
# - `h::Int64`: High orbit identifier
# """
# function getTangentBifurcationOrbit(targeter::LyapunovJCTargeter, orbitFamily::MBD.CR3BPOrbitFamily, l::Int64, h::Int64)
#     orbitl::MBD.CR3BPPeriodicOrbit = orbitFamily.familyMembers[l]
#     orbith::MBD.CR3BPPeriodicOrbit = orbitFamily.familyMembers[h]
#     newInitialCondition::Vector{Float64} = orbitl.initialCondition+0.5*(orbith.initialCondition-orbitl.initialCondition)
#     newPeriod::Float64 = orbitl.period+0.5*(orbith.period-orbitl.period)
#     newJC::Float64 = orbitl.JacobiConstant+0.5*(orbith.JacobiConstant-orbitl.JacobiConstant)
#     solutionBisect::MBD.MultipleShooterProblem = correct(targeter, newInitialCondition, [0, newPeriod], newJC, 1E-9)
#     orbitBisect = MBD.CR3BPPeriodicOrbit(solutionBisect, targeter)
#     getProperties!(targeter, orbitBisect)
#     getMonodromy!(targeter, orbitBisect)
#     getStability!(orbitBisect)
#     paramValuel::Float64 = 2*orbitl.BrouckeStability[1]+orbitl.BrouckeStability[2]+2
#     paramValueBisect::Float64 = 2*orbitBisect.BrouckeStability[1]+orbitBisect.BrouckeStability[2]+2
#     currentError = abs(paramValueBisect)
#     println("Current parameter value = $paramValueBisect")
#     counter::Int64 = 1
#     while (currentError > 1E-8) && (counter < 50)
#         (sign(paramValueBisect) == sign(paramValuel)) ? (orbitl = orbitBisect) : (orbith = orbitBisect)
#         newInitialCondition = orbitl.initialCondition+0.5*(orbith.initialCondition-orbitl.initialCondition)
#         newPeriod = orbitl.period+0.5*(orbith.period-orbitl.period)
#         newJC = orbitl.JacobiConstant+0.5*(orbith.JacobiConstant-orbitl.JacobiConstant)
#         solutionBisect = correct(targeter, newInitialCondition, [0, newPeriod], newJC, 1E-9)
#         orbitBisect = MBD.CR3BPPeriodicOrbit(solutionBisect, targeter)
#         getProperties!(targeter, orbitBisect)
#         getMonodromy!(targeter, orbitBisect)
#         getStability!(orbitBisect)
#         paramValuel = 2*orbitl.BrouckeStability[1]+orbitl.BrouckeStability[2]+2
#         paramValueBisect = 2*orbitBisect.BrouckeStability[1]+orbitBisect.BrouckeStability[2]+2
#         currentError = abs(paramValueBisect)
#         println("Current parameter value = $paramValueBisect")
#         counter += 1
#     end
#     (currentError > 1E-8) ? println("Bisection failed to converge after 50 iterations") : (return orbitBisect)
# end

"""
    propagateState(targeter, q0, tSpan)

Return propagated state

# Arguments
- `targeter::Lyapunovx0Targeter`: CR3BP Lyapunov initial x-state targeter object
- `q0::Vector{Float64}`: Initial conditions [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::Lyapunovx0Targeter, q0::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, q0, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, getStateCount(arc))
end
