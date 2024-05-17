"""
Jacobi constant perpendicular crossing targeter for CR3BP Lyapunov orbits

Author: Jonathan Richmond
C: 1/11/23
U: 5/16/24
"""

using MBD, CSV, DataFrames, LinearAlgebra

export LyapunovJCTargeter
export correct, getMonodromy!, getOrbit, getProperties!
export getTangentBifurcationOrbit, propagateState

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
    correct(targeter, q0, tSpan, targetJC; tol)

Return corrected multiple shooter problem

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Lyapunov Jacobi constant targeter object
- `q0::Vector{Float64}`: Initial conditions [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
- `targetJC::Float64`: Target Jacobi constant
- `tol::Float64`: Convergence tolerance (optional)
"""
function correct(targeter::LyapunovJCTargeter, q0::Vector{Float64}, tSpan::Vector{Float64}, targetJC::Float64, tol::Float64 = 1E-10)
    halfPeriod::Float64 = (tSpan[2]-tSpan[1])/2
    qf::Vector{Float64} = propagateState(targeter, q0, [0, halfPeriod])
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
    shooter = MBD.MultipleShooter(tol)
    #shooter.printProgress = true
    solved::MBD.MultipleShooterProblem = MBD.solve!(shooter, problem)
    convergenceCheck::MBD.ConstraintVectorL2NormConvergenceCheck = shooter.convergenceCheck
    for constraint::MBD.AbstractConstraint in getConstraints(solved)
        (LinearAlgebra.norm(evaluateConstraint(constraint, getFreeVariableIndexMap!(solved), getFreeVariableVector!(solved))) > convergenceCheck.maxVectorNorm) && println("ERROR: Multiple shooter failed to converge for $(typeof(constraint))")
    end

    return solved
end

"""
    getMonodromy!(targeter, orbit)

Return periodic orbit object with updated monodromy matrix

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Lyapunov Jacobi constant targeter object
- `orbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getMonodromy!(targeter::LyapunovJCTargeter, orbit::MBD.CR3BPPeriodicOrbit)
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    initialGuess::Vector{Float64} = appendExtraInitialConditions(targeter.dynamicsModel, orbit.problem.nodes[1].state.data, MBD.STM)
    arc::MBD.Arc = propagate(propagator, initialGuess, [0, orbit.problem.segments[1].TOF.data[1]], targeter.dynamicsModel)
    halfPeriodState::Vector{Float64} = getStateByIndex(arc, -1)
    halfPeriodSTM::Matrix{Float64} = [halfPeriodState[7:12] halfPeriodState[13:18] halfPeriodState[19:24] halfPeriodState[25:30] halfPeriodState[31:36] halfPeriodState[37:42]]
    G::Matrix{Float64} = [1.0 0 0 0 0 0; 0 -1.0 0 0 0 0; 0 0 1.0 0 0 0; 0 0 0 -1.0 0 0; 0 0 0 0 1.0 0; 0 0 0 0 0 -1.0]
    Omega::Matrix{Float64} = [0 1.0 0; -1.0 0 0; 0 0 0]
    A::Matrix{Float64} = [zeros(Float64, (3, 3)) [-1.0 0 0; 0 -1.0 0; 0 0 -1.0]; [1.0 0 0; 0 1.0 0; 0 0 1.0] -2*Omega]
    B::Matrix{Float64} = [-2*Omega [1.0 0 0; 0 1.0 0; 0 0 1.0]; [-1.0 0 0; 0 -1.0 0; 0 0 -1.0] zeros(Float64, (3, 3))]
    orbit.monodromy = G*A*(halfPeriodSTM')*B*G*halfPeriodSTM
end

"""
    getOrbit(targeter, fileName, paramName, paramValue; tol)

 Return requested periodic orbit object through bisection method

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Lyapunov Jacobi constant targeter object
- `fileName::String`: Family data file
- `paramName::String`: Requested parameter name
- `paramValue::Float64`: Requested parameter value
- `tol::Float64`: Convergence tolerance (optional)
"""
function getOrbit(targeter::LyapunovJCTargeter, filename::String, paramName::String, paramValue::Float64, tol::Float64 = 1E-10)
    familyData::DataFrames.DataFrame = DataFrames.DataFrame(CSV.File(filename))
    !any(occursin.(paramName, names(familyData))) && throw(ErrorException("Parameter not supported"))
    initialDirection::Float64 = familyData[2,paramName]-familyData[1,paramName]
    increase::Bool = (initialDirection < 0) ? false : true
    switch::Int64 = 0
    switchIndex::Vector{Int64} = [1]
    for d::Int64 in 3:size(familyData, 1)
        if (familyData[d,paramName]-familyData[d-1,paramName] > 0) && (increase == false)
            increase = true
            switch += 1
            push!(switchIndex, d-1)
        elseif (familyData[d,paramName]-familyData[d-1,paramName] < 0) && (increase == true)
            increase = false
            switch += 1
            push!(switchIndex, d-1)
        end
    end
    push!(switchIndex, size(familyData, 1))
    solutionIndex::Vector{Int64} = []
    for s::Int64 in 2:length(switchIndex)
        switchData::DataFrames.DataFrame = familyData[switchIndex[s-1]:switchIndex[s],:]
        closeIndex::Int64 = argmin(abs.(switchData[:,paramName].-paramValue))
        (abs(switchData[closeIndex,paramName]-paramValue) <= 5E-2) && push!(solutionIndex, closeIndex+switchIndex[s-1]-1)
    end
    index::Int64 = solutionIndex[1]
    if length(solutionIndex) > 1
        println("Multiple possible orbits; choosing most unstable...")
        stability::Vector{Float64} = zeros(Float64, length(solutionIndex))
        for s::Int64 in 1:length(solutionIndex)
            stability[s] = familyData[solutionIndex[s],"Stability Index"]
        end
        maxIndex::Int64 = argmax(stability)
        index = solutionIndex[maxIndex]
    end
    if abs(familyData[index,paramName]-paramValue) <= 1E-9
        println("Orbit exists in database!")
        orbitData::DataFrames.DataFrameRow = familyData[index,:]
        initialCondition::Vector{Float64} = [orbitData[s] for s in ["x", "y", "z", "xdot", "ydot", "zdot"]]
        solution = correct(targeter, initialCondition, [0, orbitData["Period"]], orbitData["JC"], tol)
        orbitBisect = MBD.CR3BPPeriodicOrbit(solution, targeter)
        getProperties!(targeter, orbitBisect)
        getMonodromy!(targeter, orbitBisect)
        getStability!(orbitBisect)
        return orbitBisect
    else
        println("Orbit does not already exist in database; attempting bisection...")
        subData::DataFrames.DataFrame = DataFrames.sort(familyData[index-1:index+1,:], paramName)
        index2::Int64 = 3
        ((subData[2,paramName] - paramValue) > 0) && (index2 = 1)
        data1::DataFrames.DataFrameRow = subData[2,:]
        data2::DataFrames.DataFrameRow = subData[index2,:]
        newInitialCondition::Vector{Float64} = [data1[s]+0.5*(data2[s]-data1[s]) for s = ["x", "y", "z", "xdot", "ydot", "zdot"]]
        newPeriod::Float64 = data1["Period"]+0.5*(data2["Period"]-data1["Period"])
        newJC::Float64 = data1["JC"]+0.5*(data2["JC"]-data1["JC"])
        currentError::Float64 = abs(data1[paramName]-paramValue)
        counter::Int64 = 0
        while (currentError > 1E-9) && (counter < 50)
            solutionBisect::MBD.MultipleShooterProblem = correct(targeter, newInitialCondition, [0, newPeriod], newJC, tol)
            orbitBisect = MBD.CR3BPPeriodicOrbit(solutionBisect, targeter)
            getProperties!(targeter, orbitBisect)
            getMonodromy!(targeter, orbitBisect)
            getStability!(orbitBisect)
            dataBisect::DataFrames.DataFrameRow = DataFrames.DataFrame("x" => orbitBisect.initialCondition[1], "y" => orbitBisect.initialCondition[2], "z" => orbitBisect.initialCondition[3], "xdot" => orbitBisect.initialCondition[4], "ydot" => orbitBisect.initialCondition[5], "zdot" => orbitBisect.initialCondition[6], "JC" => orbitBisect.JacobiConstant, "Period" => orbitBisect.period, "Stability Index" => orbitBisect.nu)[1,:]
            currentError = abs(dataBisect[paramName]-paramValue)
            println("Current parameter value: "*paramName*" = $(dataBisect[paramName])")
            (sign(dataBisect[paramName]-paramValue) == sign(data1[paramName]-paramValue)) ? (data1 = dataBisect) : (data2 = dataBisect)
            newInitialCondition = [data1[s]+0.5*(data2[s]-data1[s]) for s in ["x", "y", "z", "xdot", "ydot", "zdot"]]
            newPeriod = data1["Period"]+0.5*(data2["Period"]-data1["Period"])
            newJC = data1["JC"]+0.5*(data2["JC"]-data1["JC"])
            counter += 1
        end
        (currentError > 1E-9) ? println("Bisection failed to converge after 50 iterations") : println("Bisection successful!")
        return orbitBisect
    end
end

"""
    getProperties!(targeter, orbit)

Return periodic orbit object with updated properties

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Lyapunov Jacobi constant targeter object
- `orbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getProperties!(targeter::LyapunovJCTargeter, orbit::MBD.CR3BPPeriodicOrbit)
    initialState::Vector{Float64} = orbit.problem.nodes[1].state.data
    period::Float64 = 2*orbit.problem.segments[1].TOF.data[1]
    finalState::Vector{Float64} = propagateState(targeter, initialState, [0.0, period])
    (LinearAlgebra.norm(finalState-initialState) > 1E-6) && throw(ErrorException("Trajectory not periodic (error norm = $(LinearAlgebra.norm(finalState-initialState)))"))
    orbit.initialCondition = initialState
    orbit.period = period
    orbit.JacobiConstant = getJacobiConstant(targeter.dynamicsModel, initialState)
end

"""
    getTangentBifurcationOrbit(targeter, orbitFamily, l, h)

Return tangent bifurcation orbit

# Arguments
- `targeter::LyapunovJCTargeter`: CR3BP Laypunov Jacobi constant targeter object
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
- `l::Int64`: Low orbit identifier
- `h::Int64`: High orbit identifier
"""
function getTangentBifurcationOrbit(targeter::LyapunovJCTargeter, orbitFamily::MBD.CR3BPOrbitFamily, l::Int64, h::Int64)
    orbitl::MBD.CR3BPPeriodicOrbit = orbitFamily.familyMembers[l]
    orbith::MBD.CR3BPPeriodicOrbit = orbitFamily.familyMembers[h]
    newInitialCondition::Vector{Float64} = orbitl.initialCondition+0.5*(orbith.initialCondition-orbitl.initialCondition)
    newPeriod::Float64 = orbitl.period+0.5*(orbith.period-orbitl.period)
    newJC::Float64 = orbitl.JacobiConstant+0.5*(orbith.JacobiConstant-orbitl.JacobiConstant)
    solutionBisect::MBD.MultipleShooterProblem = correct(targeter, newInitialCondition, [0, newPeriod], newJC, 1E-9)
    orbitBisect = MBD.CR3BPPeriodicOrbit(solutionBisect, targeter)
    getProperties!(targeter, orbitBisect)
    getMonodromy!(targeter, orbitBisect)
    getStability!(orbitBisect)
    paramValuel::Float64 = 2*orbitl.BrouckeStability[1]+orbitl.BrouckeStability[2]+2
    paramValueBisect::Float64 = 2*orbitBisect.BrouckeStability[1]+orbitBisect.BrouckeStability[2]+2
    currentError = abs(paramValueBisect)
    println("Current parameter value = $paramValueBisect")
    counter::Int64 = 1
    while (currentError > 1E-8) && (counter < 50)
        (sign(paramValueBisect) == sign(paramValuel)) ? (orbitl = orbitBisect) : (orbith = orbitBisect)
        newInitialCondition = orbitl.initialCondition+0.5*(orbith.initialCondition-orbitl.initialCondition)
        newPeriod = orbitl.period+0.5*(orbith.period-orbitl.period)
        newJC = orbitl.JacobiConstant+0.5*(orbith.JacobiConstant-orbitl.JacobiConstant)
        solutionBisect = correct(targeter, newInitialCondition, [0, newPeriod], newJC, 1E-9)
        orbitBisect = MBD.CR3BPPeriodicOrbit(solutionBisect, targeter)
        getProperties!(targeter, orbitBisect)
        getMonodromy!(targeter, orbitBisect)
        getStability!(orbitBisect)
        paramValuel = 2*orbitl.BrouckeStability[1]+orbitl.BrouckeStability[2]+2
        paramValueBisect = 2*orbitBisect.BrouckeStability[1]+orbitBisect.BrouckeStability[2]+2
        currentError = abs(paramValueBisect)
        println("Current parameter value = $paramValueBisect")
        counter += 1
    end
    (currentError > 1E-8) ? println("Bisection failed to converge after 50 iterations") : (return orbitBisect)
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
