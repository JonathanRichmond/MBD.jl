"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 9/1/22
U: 5/15/24
"""

using MBD
using Test

include("../src/bifurcation/Bifurcation.jl")
include("../src/continuation/AdaptiveStepSizeByElementGenerator.jl")
include("../src/continuation/BoundingBoxContinuationEndCheck.jl")
include("../src/continuation/BoundingBoxJumpCheck.jl")
include("../src/continuation/JacobiConstantContinuationEngine.jl")
include("../src/continuation/NaturalParameterContinuationEngine.jl")
include("../src/continuation/NumberStepsContinuationEndCheck.jl")
include("../src/corrections/ConstraintVectorL2NormConvergenceCheck.jl")
include("../src/corrections/ContinuityConstraint.jl")
include("../src/corrections/LeastSquaresUpdateGenerator.jl")
include("../src/corrections/MinimumNormUpdateGenerator.jl")
include("../src/corrections/MultipleShooter.jl")
include("../src/corrections/MultipleShooterProblem.jl")
include("../src/corrections/Node.jl")
include("../src/corrections/Segment.jl")
include("../src/corrections/StateConstraint.jl")
include("../src/corrections/StateMatchConstraint.jl")
include("../src/corrections/Variable.jl")
include("../src/CR3BP/DynamicsModel.jl")
include("../src/CR3BP/EquationsOfMotion.jl")
include("../src/CR3BP/JacobiConstraint.jl")
include("../src/CR3BP/OrbitFamily.jl")
include("../src/CR3BP/PeriodicOrbit.jl")
include("../src/CR3BP/SystemData.jl")
include("../src/propagation/Arc.jl")
include("../src/propagation/EventFunctions.jl")
include("../src/propagation/Propagator.jl")
include("../src/spice/BodyName.jl")
include("../src/spice/SpiceFunctions.jl")
include("../src/TBP/DynamicsModel.jl")
include("../src/TBP/EquationsOfMotion.jl")
include("../src/TBP/Trajectory.jl")
include("../src/utilities/UtilityFunctions.jl")

include("ExampleLyapunovJCTargeter.jl")

@testset "Files" begin
    @test isfile("../src/body_data.xml")
    @test isfile("../src/spice/kernels/de440.bsp")
    @test isfile("../src/spice/kernels/naif0012.tls")
end

@testset "Constructors" begin
    Earth = MBD.BodyName("Earth")
    @test Earth.name == "Earth"
    EarthData = MBD.BodyData("Earth")
    @test EarthData.spiceID == 399
    @test EarthData.mass == 5.972580035423323E24
    @test MBD.CR3BPSystemData <: MBD.AbstractSystemData
    @test_throws ArgumentError MBD.CR3BPSystemData("Earth", "Mars")
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    @test systemData.primarySpiceIDs == [399, 301]
    @test systemData.params[1] == 0.012150584269940356
    @test MBD.CR3BPDynamicsModel <: MBD.AbstractDynamicsModel
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    @test dynamicsModel.systemData == systemData
    @test MBD.CR3BPEquationsOfMotion <: MBD.AbstractEquationsOfMotion
    EOMs = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    @test EOMs.dim == 6
    @test EOMs.mu == 0.012150584269940356
    integratorFactory = MBD.IntegratorFactory(MBD.DP8)
    @test integratorFactory.integrator == DifferentialEquations.DP8()
    propagator = MBD.Propagator()
    @test propagator.absTol == 1E-12
    propagatorBS = MBD.Propagator(MBD.BS5)
    @test propagatorBS.integratorFactory == MBD.IntegratorFactory(MBD.BS5)
    arc = MBD.Arc(dynamicsModel)
    @test arc.dynamicsModel == dynamicsModel
    @test_throws ArgumentError MBD.Variable([1.0, 1.0], [true, false, true])
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test stateVariable.freeVariableMask == [true, false, false, false, true, false]
    @test MBD.Node <: MBD.IHasVariables
    @test_throws ArgumentError MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263], dynamicsModel)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    @test originNode.state.data == stateVariable.data
    @test originNode.state.name == "Node State"
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    @test MBD.Segment <: MBD.IHasVariables
    @test_throws ArgumentError MBD.Segment(2.743, originNode, originNode)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    @test segment.TOF == MBD.Variable([2.743], [true])
    @test MBD.MultipleShooterProblem <: MBD.AbstractNonlinearProblem
    problem = MBD.MultipleShooterProblem()
    @test !problem.hasBeenBuilt
    @test MBD.ContinuityConstraint <: MBD.AbstractConstraint
    continuityConstraint = MBD.ContinuityConstraint(segment)
    @test continuityConstraint.constrainedIndices == 1:6
    @test MBD.StateConstraint <: MBD.AbstractConstraint
    @test_throws ArgumentError MBD.StateConstraint(originNode, [1, 5], [1.0])
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    @test_throws ArgumentError MBD.StateConstraint(originNode, [2], [1.0])
    stateConstraint = MBD.StateConstraint(originNode, [1, 5], [0.8, 0.1])
    @test stateConstraint.variable == originNode.state
    @test MBD.StateMatchConstraint <: MBD.AbstractConstraint
    variable1 = MBD.Variable([1.0, 2.0], [true, false])
    setFreeVariableMask!(terminalNode.state, [true, false, true, true, true, true])
    @test_throws ArgumentError MBD.StateMatchConstraint(stateVariable, variable1, [1, 2])
    @test_throws ArgumentError MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2])
    setFreeVariableMask!(terminalNode.state, [true, true, true, true, true, true])
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test stateMatchConstraint.variable1 == originNode.state
    @test MBD.JacobiConstraint <: MBD.AbstractConstraint
    systemDataTBP = MBD.TBPSystemData("Earth")
    dynamicsModelTBP = MBD.TBPDynamicsModel(systemDataTBP)
    nodeTBP = MBD.Node(0.0, originNode.state.data, dynamicsModelTBP)
    @test_throws ArgumentError MBD.JacobiConstraint(nodeTBP, 3.0)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    @test JCConstraint.epoch == originNode.epoch
    @test MBD.ConstraintVectorL2NormConvergenceCheck <: MBD.AbstractConvergenceCheck
    convergenceCheck = MBD.ConstraintVectorL2NormConvergenceCheck()
    @test convergenceCheck.maxVectorNorm == 1E-10
    convergenceCheckTight = MBD.ConstraintVectorL2NormConvergenceCheck(1E-12)
    @test convergenceCheckTight.maxVectorNorm == 1E-12
    @test MBD.MinimumNormUpdateGenerator <: MBD.AbstractUpdateGenerator
    @test MBD.LeastSquaresUpdateGenerator <: MBD.AbstractUpdateGenerator
    @test MBD.MultipleShooter <: MBD.AbstractNonlinearProblemSolver
    multipleShooter = MBD.MultipleShooter()
    @test multipleShooter.maxIterations == 25
    multipleShooterTight = MBD.MultipleShooter(1E-12)
    @test multipleShooterTight.convergenceCheck.maxVectorNorm == 1E-12
    continuationData = MBD.ContinuationData()
    @test continuationData.converging
    @test_throws ArgumentError MBD.AdaptiveStepSizeByElementGenerator("Initial State", -1, 1E-4, 1E-2)
    stepSizeGenerator = MBD.AdaptiveStepSizeByElementGenerator("Initial State", 1, -1E-4, -1E-2)
    @test stepSizeGenerator.minStepSize == -1E-10
    @test MBD.NumberStepsContinuationEndCheck <: MBD.AbstractContinuationEndCheck
    stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
    @test stepsEndCheck.maxSteps == 100
    @test MBD.BoundingBoxContinuationEndCheck <: MBD.AbstractContinuationEndCheck
    boxEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [0 0.99; NaN NaN])
    @test isequal(boxEndCheck.paramBounds, [0 0.99; NaN NaN])
    @test MBD.BoundingBoxJumpCheck <: MBD.AbstractContinuationJumpCheck
    boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
    @test isequal(boxJumpCheck.paramBounds, [0 1.2; 0 1.0])
    @test MBD.NaturalParameterContinuationEngine <: MBD.AbstractContinuationEngine
    naturalParamContinuation = MBD.NaturalParameterContinuationEngine("Initial State", 1, -1E-4, -1E-2)
    @test naturalParamContinuation.storeIntermediateMembers
    naturalParamContinuationTight = MBD.NaturalParameterContinuationEngine("Initial State", 1, -1E-4, -1E-2, 1E-12)
    @test naturalParamContinuationTight.corrector.convergenceCheck.maxVectorNorm == 1E-12
    @test MBD.JacobiConstantContinuationEngine <: MBD.AbstractContinuationEngine
    JCContinuation = MBD.JacobiConstantContinuationEngine(-1E-4, -1E-1)
    @test JCContinuation.storeIntermediateMembers
    JCContinuationTight = MBD.JacobiConstantContinuationEngine(-1E-4, -1E-1, 1E-12)
    @test JCContinuationTight.corrector.convergenceCheck.maxVectorNorm == 1E-12
    #bifurcation = MBD.Bifurcation(orbitFamily, periodicOrbit, 1, MBD.TANGENT, 1)
    #@test bifurcation.orbit == 1
    @test MBD.CR3BPPeriodicOrbit <: MBD.AbstractTrajectoryStructure
    targeter = LyapunovJCTargeter(dynamicsModel)
    periodicOrbit = MBD.CR3BPPeriodicOrbit(problem, targeter)
    @test periodicOrbit.tau == 0
    @test MBD.CR3BPOrbitFamily <: MBD.AbstractStructureFamily
    orbitFamily = MBD.CR3BPOrbitFamily([periodicOrbit])
    @test length(orbitFamily.familyMembers) == 1
    @test MBD.CR3BPManifoldArc <: MBD.AbstractTrajectoryStructure
    manifoldArc = MBD.CR3BPManifoldArc([0.8234+0im, 0, 0, 0, 0.1263+0im, 0], periodicOrbit)
    @test manifoldArc.TOF == 0
    @test MBD.TBPSystemData <: MBD.AbstractSystemData
    @test systemDataTBP.charTime == 2.898139774118648e9
    @test MBD.TBPDynamicsModel <: MBD.AbstractDynamicsModel
    @test dynamicsModelTBP.systemData == systemDataTBP
    @test MBD.TBPEquationsOfMotion <: MBD.AbstractEquationsOfMotion
    EOMsTBP = MBD.TBPEquationsOfMotion(MBD.SIMPLE, dynamicsModelTBP)
    @test EOMsTBP.mu == systemDataTBP.gravParam
    @test MBD.TBPTrajectory <: MBD.AbstractTrajectoryStructure
    trajectory = MBD.TBPTrajectory(originNode.state.data, dynamicsModelTBP)
    @test trajectory.E == 0
end

@testset "Copy" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    @test MBD.shallowClone(originNode) == originNode
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    @test MBD.shallowClone(segment) == segment
    problem = MBD.MultipleShooterProblem()
    @test MBD.shallowClone(problem) == problem
    continuityConstraint = MBD.ContinuityConstraint(segment)
    @test MBD.shallowClone(continuityConstraint) == continuityConstraint
    stateConstraint = MBD.StateConstraint(originNode, [1, 5], [0.8, 0.1])
    @test MBD.shallowClone(stateConstraint) == stateConstraint
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test MBD.shallowClone(stateMatchConstraint) == stateMatchConstraint
    jacobiConstraint = MBD.JacobiConstraint(originNode, 3.0)
    @test MBD.shallowClone(jacobiConstraint) == jacobiConstraint
    targeter = LyapunovJCTargeter(dynamicsModel)
    periodicOrbit = MBD.CR3BPPeriodicOrbit(problem, targeter)
    @test shallowClone(periodicOrbit) == periodicOrbit
    systemDataTBP = MBD.TBPSystemData("Earth")
    dynamicsModelTBP = MBD.TBPDynamicsModel(systemDataTBP)
    trajectory = MBD.TBPTrajectory(originNode.state.data, dynamicsModelTBP)
    @test shallowClone(trajectory) == trajectory
end

@testset "Deep Copy" begin
    problem = MBD.MultipleShooterProblem()
    @test MBD.deepClone(problem) == problem
    variable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test MBD.deepClone(variable) == variable
end

#@testset "Bifurcation" begin    
#end

@testset "BodyName" begin
    bodyName = MBD.BodyName("Earth")
    @test getIDCode(bodyName) == 399
end

@testset "AdaptiveStepSizeByElementGenerator" begin
    stepSizeGenerator = MBD.AdaptiveStepSizeByElementGenerator("Initial State", 1, -1E-4, -1E-2)
    continuationData = MBD.ContinuationData()
    continuationData.currentStepSize = stepSizeGenerator.initialStepSize
    updateStepSize!(stepSizeGenerator, continuationData)
    @test continuationData.currentStepSize == -0.0002
    continuationData.numIterations = 15
    updateStepSize!(stepSizeGenerator, continuationData)
    @test continuationData.currentStepSize == -0.0001
    continuationData.numIterations = 1
    continuationData.currentStepSize = -0.05
    updateStepSize!(stepSizeGenerator, continuationData)
    @test continuationData.currentStepSize == -0.01
    continuationData.currentStepSize = stepSizeGenerator.initialStepSize
    continuationData.converging = false
    updateStepSize!(stepSizeGenerator, continuationData)
    @test continuationData.currentStepSize == -0.00005
    continuationData.currentStepSize = -1E-10
    updateStepSize!(stepSizeGenerator, continuationData)
    @test continuationData.forceEndContinuation
end

@testset "BoundingBoxContinuationEndCheck" begin
    boxEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [0 0.99; NaN NaN])
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    stateVariable.name = "Initial State"
    checkBounds(boxEndCheck, stateVariable)
    setFreeVariableMask!(stateVariable, [true, false, false, false, false, false])
    @test_throws ArgumentError checkBounds(boxEndCheck, stateVariable)
    boxEndCheck.paramBounds = [0 0.5 0.99]
    @test_throws ArgumentError checkBounds(boxEndCheck, stateVariable)
    boxEndCheck.paramBounds = [0.5 0]
    @test_throws ArgumentError checkBounds(boxEndCheck, stateVariable)
    #isContinuationDone
end

@testset "BoundingBoxJumpCheck" begin
    boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 0.99; NaN NaN])
    problem = MBD.MultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    addBounds!(boxJumpCheck, problem, problem.segments[1].originNode.state, [0 0.9; NaN NaN])
    @test length(boxJumpCheck.variableBounds) == 1
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test_throws MethodError addBounds!(boxJumpCheck, problem, stateVariable, [0 0.89; NaN NaN])
    setFreeVariableMask!(stateVariable, [true, false, false, false, false, false])
    @test_throws ArgumentError checkBounds(boxJumpCheck, stateVariable, [0 0.9; NaN NaN])
    @test_throws ArgumentError checkBounds(boxJumpCheck, stateVariable, [0 0.5 0.9])
    @test_throws ArgumentError checkBounds(boxJumpCheck, stateVariable, [0.9 0])
    removeBounds!(boxJumpCheck, problem, problem.segments[1].originNode.state)
    @test length(boxJumpCheck.variableBounds) == 0
    @test_throws MethodError removeBounds!(boxJumpCheck, problem, stateVariable)
    #isFamilyMember
end

@testset "JacobiConstantContinuationEngine" begin
    JCContinuation = MBD.JacobiConstantContinuationEngine(-1E-4, -1E-1)
    stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
    addEndCheck!(JCContinuation, stepsEndCheck)
    @test length(JCContinuation.endChecks) == 1
    boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
    addJumpCheck!(JCContinuation, boxJumpCheck)
    @test length(JCContinuation.jumpChecks) == 1
    JCContinuation.dataInProgress.converging = false
    resetEngine!(JCContinuation)
    @test JCContinuation.dataInProgress.converging
    #computeFullStep
    #constrainNextGuess!
    #convergeInitialSolution
    #doContinuation!
    #endContinuation
    #tryConverging!
end

@testset "NaturalParameterContinuationEngine" begin
    naturalParamContinuation = MBD.NaturalParameterContinuationEngine("Initial State", 1, -1E-4, -1E-2)
    stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
    addEndCheck!(naturalParamContinuation, stepsEndCheck)
    @test length(naturalParamContinuation.endChecks) == 1
    boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
    addJumpCheck!(naturalParamContinuation, boxJumpCheck)
    @test length(naturalParamContinuation.jumpChecks) == 1
    naturalParamContinuation.dataInProgress.converging = false
    resetEngine!(naturalParamContinuation)
    @test naturalParamContinuation.dataInProgress.converging
    #computeFullStep
    #constrainNextGuess!
    #convergeInitialSolution
    #doContinuation!
    #endContinuation
    #tryConverging!
end

@testset "NumberStepsContinuationEndCheck" begin
    stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
    continuationData = MBD.ContinuationData()
    @test !isContinuationDone(stepsEndCheck, continuationData)
    continuationData.stepCount = 100
    oldstd = stdout
    redirect_stdout(devnull)
    @test isContinuationDone(stepsEndCheck, continuationData)
    redirect_stdout(oldstd)
end

@testset "ConstraintVectorL2NormConvergenceCheck" begin
    convergenceCheck = MBD.ConstraintVectorL2NormConvergenceCheck()
    problem = MBD.MultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.Node(2.743, [0.8, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    addConstraint!(problem, continuityConstraint)
    @test !isConverged(convergenceCheck, problem)
end

@testset "ContinuityConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    problem = MBD.MultipleShooterProblem()
    @test isApproxSigFigs(evaluateConstraint(continuityConstraint, problem.freeVariableIndexMap, problem.freeVariableVector), [-1.535771509964E-11, 4.571839868489E-12, 0, -4.366933897826E-11, 1.664283988401E-11, 0], 13)
    @test getNumberConstraintRows(continuityConstraint) == 6
    @test length(keys(getPartials_ConstraintWRTVariables(continuityConstraint, problem.freeVariableIndexMap, problem.freeVariableVector))) == 3
    #updatePointers!
end

@testset "LeastSquaresUpdateGenerator" begin
    leastSquaresUpdate = MBD.LeastSquaresUpdateGenerator()
    problem = MBD.MultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    addConstraint!(problem, continuityConstraint)
    @test !canGenerateUpdate(leastSquaresUpdate, problem)
    @test_throws ErrorException getFullUpdate(leastSquaresUpdate, problem)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    stateConstraint = MBD.StateConstraint(terminalNode, [2, 4, 6], [0.0, 0.0, 0.0])
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, stateConstraint)
    @test canGenerateUpdate(leastSquaresUpdate, problem)
    @test isApproxSigFigs([LinearAlgebra.norm(getFullUpdate(leastSquaresUpdate, problem))], [1.2820039], 8)
end

@testset "MinimumNormUpdateGenerator" begin
    minimumNormUpdate = MBD.MinimumNormUpdateGenerator()
    problem = MBD.MultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    stateConstraint = MBD.StateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(problem, continuityConstraint)
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, stateConstraint)
    @test canGenerateUpdate(minimumNormUpdate, problem)
    @test isApproxSigFigs([LinearAlgebra.norm(getFullUpdate(minimumNormUpdate, problem))], [1.2820039], 8)
    initialStateConstraint = MBD.StateConstraint(originNode, [1], [0.82])
    addConstraint!(problem, initialStateConstraint)
    @test !canGenerateUpdate(minimumNormUpdate, problem)
    @test_throws ErrorException getFullUpdate(minimumNormUpdate, problem)
end

@testset "MultipleShooter" begin
    multipleShooter = MBD.MultipleShooter()
    setPrintProgress!(multipleShooter, true)
    @test multipleShooter.printProgress
    setPrintProgress!(multipleShooter, false)
    problem = MBD.MultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    stateConstraint = MBD.StateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(problem, continuityConstraint)
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, stateConstraint)
    @test_throws ErrorException solve!(multipleShooter, problem)
    newProblem = MBD.MultipleShooterProblem()
    addSegment!(newProblem, segment)
    newJCConstraint = MBD.JacobiConstraint(originNode, 3.1743)
    addConstraint!(newProblem, continuityConstraint)
    addConstraint!(newProblem, newJCConstraint)
    addConstraint!(newProblem, stateConstraint)
    @test isApproxSigFigs([LinearAlgebra.norm(solveUpdateEquation(multipleShooter, newProblem))], [0.027973617], 8)
    solution::MBD.MultipleShooterProblem = solve!(multipleShooter, newProblem)
    @test isConverged(multipleShooter.convergenceCheck, solution)
end

@testset "MultipleShooterProblem" begin
    problem = MBD.MultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    @test length(problem.segments) == 1
    @test_throws ErrorException buildAdjacencyMatrix(problem)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    addConstraint!(problem, continuityConstraint)
    @test length(problem.constraintIndexMap) == 1
    buildProblem!(problem)
    @test problem.hasBeenBuilt
    @test checkJacobian(problem)
    @test length(getConstraints(problem)) == 1
    @test isApproxSigFigs([LinearAlgebra.norm(getConstraintVector!(problem))], [4.804712395091E-5], 13)
    @test length(getFreeVariableIndexMap!(problem)) == 6
    @test isApproxSigFigs([LinearAlgebra.norm(getFreeVariableVector!(problem))], [2.987433787532], 13)
    @test size(getJacobian(problem)) == (6, 9)
    @test getNumberConstraints(problem) == 6
    @test getNumberFreeVariables(problem) == 9
    removeConstraint!(problem, collect(keys(problem.constraintIndexMap))[1])
    @test length(problem.constraintIndexMap) == 0
    newFreeVariableVector::Vector{Float64} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    setFreeVariableVector!(problem, newFreeVariableVector)
    @test getFreeVariableVector!(problem) == newFreeVariableVector
    badProblem = MBD.MultipleShooterProblem()
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    addVariable!(badProblem, stateVariable)
    @test length(badProblem.freeVariableIndexMap) == 1
    #checkValidGraph
end

@testset "Node" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    node = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, true, true, true, true, true])
    stateVariable.name = "Node State"
    epochVariable = MBD.Variable([0.0], [false])
    epochVariable.name = "Node Epoch"
    @test getVariables(node) == [stateVariable, epochVariable]
    #updatePointers!
end

@testset "Segment" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    @test isApproxSigFigs(getFinalState!(segment), [0.8322038365943, -0.002795978474764, 0, 0.02460934345210, 0.1166002944939, 0], 13)
    @test isApproxSigFigs(getFinalStateRate!(segment), [0.02460934345210, 0.1166002944939, 0, 0.1811347722610, -0.03842090358821, 0], 13)
    @test getPartials_FinalStateWRTEpoch!(segment) == zeros(Float64, (6, 1))
    @test isApproxSigFigs(getPartials_FinalStateWRTInitialState!(segment)[1,:], [1363.762346338, -350.5189329933, 0, 401.4700668419, 130.9118799090, 0], 13)
    @test size(getPartials_FinalStateWRTParams!(segment)) == (6, 0)
    @test getPropagatorParametersData(segment) == []
    TOFVariable = MBD.Variable([2.743], [true])
    TOFVariable.name = "Segment Time-of-Flight"
    paramsVariable = MBD.Variable(Array{Float64}(undef, 0), Array{Bool}(undef, 0))
    paramsVariable.name = "Segment Propagation Parameters"
    @test getVariables(segment) == [TOFVariable, paramsVariable]
    lazyPropagate!(segment, MBD.FULL)
    @test isApproxSigFigs(getStateByIndex(segment.propArc, -1), [0.8322038366131, -0.002795978480393, 0, 0.02460934350549, 0.1166002944736, 0, 1363.762346338, -421.0811592079, 0, 3877.981886745, -1465.214687088, 0, -350.5189329933, 109.0762551667, 0, -995.8827212621, 376.8074206179, 0, 0, 0, 0.9853225289758, 0, 0, -0.09704680344079, 401.4700668419, -123.8052684905, 0, 1141.651390589, -431.6882797629, 0, 130.9118799090, -40.55137319719, 0, 372.4229726018, -139.7310684168, 0, 0, 0, -0.07418084339284, 0, 0, 1.022202359236], 13)
    resetPropagatedArc!(segment)
    println(segment.propArc.states)
    @test getStateCount(segment.propArc) == 1
    #updatePointers!
end

@testset "StateConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    stateConstraint = MBD.StateConstraint(originNode, [1], [originNode.state.data[1]])
    problem = MBD.MultipleShooterProblem()
    @test evaluateConstraint(stateConstraint, problem.freeVariableIndexMap, problem.freeVariableVector) == [0]
    @test getNumberConstraintRows(stateConstraint) == 1
    @test length(keys(getPartials_ConstraintWRTVariables(stateConstraint, problem.freeVariableIndexMap, problem.freeVariableVector))) == 1
    #updatePointers!
end

@testset "StateMatchConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    problem = MBD.MultipleShooterProblem()
    @test isApproxSigFigs(evaluateConstraint(stateMatchConstraint, problem.freeVariableIndexMap, problem.freeVariableVector), [-0.008803836609691, 0.002795978479336, 0, -0.02460934349576, 0.009699705522694, 0], 13)
    @test getNumberConstraintRows(stateMatchConstraint) == 6
    @test length(keys(getPartials_ConstraintWRTVariables(stateMatchConstraint, problem.freeVariableIndexMap, problem.freeVariableVector))) == 2
end

@testset "Variable" begin
    variable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test getData(variable) == [0.8234, 0, 0, 0, 0.1263, 0]
    @test getFreeVariableData(variable) == [0.8234, 0.1263]
    @test getFreeVariableMask(variable) == [true, false, false, false, true, false]
    @test getNumberFreeVariables(variable) == 2
    @test_throws ArgumentError setFreeVariableData!(variable, [4.0, 7.0, 8.5])
    setFreeVariableData!(variable, [4.0, 7.0])
    @test variable.data == [4.0, 0, 0, 0, 7.0, 0]
    @test_throws ArgumentError setFreeVariableMask!(variable, [true, false])
    setFreeVariableMask!(variable, [true, false, false, false, false, false])
    @test variable.freeVariableMask == [true, false, false, false, false, false]
end

@testset "CR3BPSystemData" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    @test getMassRatio(systemData) == 0.012150584269940356
end

@testset "CR3BPDynamicsModel" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    @test getStateSize(dynamicsModel, MBD.SIMPLE) == 6
    @test_throws ArgumentError appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1], MBD.SIMPLE)
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.SIMPLE) == [0.8234, 0, 0, 0, 0.1263, 0]
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL) == [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    @test evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.8234, 0, 0, 0, 0.1263, 0], systemData.params) == [0, 0.1263, 0, 0.11033238649399063, 0, 0]
    MoonData = MBD.BodyData("Moon")
    @test get2BApproximation(dynamicsModel, MoonData, 2, 0.01) == [0.9778494157300597, 0, 0, 0, 1.1122968869565202, 0]
    @test getEquationsOfMotion(dynamicsModel, MBD.SIMPLE) == MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    @test isEpochIndependent(dynamicsModel)
    @test_throws ArgumentError getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [0, 0, 0, 0, 0, 0]
    @test_throws ArgumentError getEquilibriumPoint(dynamicsModel, 6)
    @test getEquilibriumPoint(dynamicsModel, 1) == [0.8369151323643023, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, 2) == [1.1556821602923404, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, 3) == [-1.005062645252109, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, 4) == [0.48784941573005963, 0.8660254037844386, 0]
    @test getJacobiConstant(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0]) == 3.1743560232059265
    L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
    L5::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 5)
    @test_throws ArgumentError getLinearVariation(dynamicsModel, 6, L1, [0.005, 0, 0])
    @test getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0]) == ([0.8419151323643023, 0, 0, 0, -0.04186136597442648, 0], [0, 2.6915795607981865])
    @test_throws ArgumentError getLinearVariation(dynamicsModel, 5, L5, [0.005, 0, 0], "Medium")
    @test getLinearVariation(dynamicsModel, 5, L5, [0.005, 0, 0], "Long") == ([0.49284941573005964, -0.8660254037844386, 0, -0.003168674904327207, -0.0020973202593643175, 0], [0, 21.06979705455942])
    @test_throws ArgumentError getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test size(getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])) == (6, 0)
    @test_throws ArgumentError getPrimaryPosition(dynamicsModel, 3)
    @test getPrimaryPosition(dynamicsModel, 1) == [-0.012150584269940356, 0, 0]
    @test getPseudopotentialJacobian(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0]) == [9.851145859594304, -3.4255729297971507, -4.42557292979715, 0, 0, 0]
    @test_throws ArgumentError getStateTransitionMatrix(dynamicsModel, [0.8324, 0, 0, 0, 0.1263, 0])
    @test getStateTransitionMatrix(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
    @test_throws ArgumentError primaryInertial2Rotating(dynamicsModel, 2, [[1.0, 0, 0, 0, 1.0, 0], [1.1, 0, 0, 0, 1.1, 0]], [0.0])
    @test_throws ArgumentError primaryInertial2Rotating(dynamicsModel, 3, [[1.0, 0, 0, 0, 1.0, 0]], [0.0])
    @test primaryInertial2Rotating(dynamicsModel, 1, [[0.8234, 0, 0, 0, 0.1263, 0]], [0.0]) == [[0.8112494157300597, 0, 0, 0, -0.6971, 0]]
    @test rotating2PrimaryInertial(dynamicsModel, 1, [[0.8234, 0, 0, 0, 0.1263, 0]], [0.0]) == [[0.8355505842699403, 0, 0, 0, 0.9618505842699403, 0]]
    SMSystemData = MBD.CR3BPSystemData("Sun", "Mars")
    SMDynamicsModel = MBD.CR3BPDynamicsModel(SMSystemData)
end

@testset "CR3BPEquationsOfMotion" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    EOMs_simple = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    qdot_simple = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.SIMPLE))
    computeDerivatives!(qdot_simple, [0.8234, 0, 0, 0, 0.1263, 0], (EOMs_simple,), 0.0)
    @test qdot_simple == [0, 0.1263, 0, 0.11033238649399063, 0, 0]
    EOMs_full = MBD.CR3BPEquationsOfMotion(MBD.FULL, dynamicsModel)
    qdot_full = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.FULL))
    computeDerivatives!(qdot_full, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL), (EOMs_full,), 0.0)
    @test qdot_full == [0, 0.1263, 0, 0.11033238649399063, 0, 0, 0, 0, 0, 9.851145859594295, 0, 0, 0, 0, 0, 0, -3.425572929797148, 0, 0, 0, 0, 0, 0, -4.4255729297971484, 1, 0, 0, 0, -2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0]
end

@testset "Propagator" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.Arc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test arc.states[end] == [0.8322038365943337, -0.0027959784747644103, 0, 0.024609343452095516, 0.11660029449394844, 0]
    @test arc.times[end-1] == 2.6859388949346146
end

@testset "Arc" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.Arc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test getStateByIndex(arc, -1) == [0.8322038365943337, -0.0027959784747644103, 0, 0.024609343452095516, 0.11660029449394844, 0]
    @test getTimeByIndex(arc, -1) == 2.743
    deleteStateAndTime!(arc, 32)
    @test getStateCount(arc) == 32
    @test arc.states[32] == [0.8322038365943337, -0.0027959784747644103, 0, 0.024609343452095516, 0.11660029449394844, 0]
    @test arc.times[32] == 2.743
    @test_throws BoundsError deleteStateAndTime!(arc, 33)
    @test_throws BoundsError getStateByIndex(arc, 33)
    @test_throws BoundsError getTimeByIndex(arc, 33)
    setParameters!(arc, systemData.params)
    @test arc.params == systemData.params
end

@testset "JacobiConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    jacobiConstraint = MBD.JacobiConstraint(originNode, 3.1743560232059265)
    @test getNumberConstraintRows(jacobiConstraint) == 1
    multipleShooterProblem = MBD.MultipleShooterProblem()
    @test evaluateConstraint(jacobiConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector) == [0.0]
    @test length(keys(getPartials_ConstraintWRTVariables(jacobiConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector))) == 1
end

@testset "TBPDynamicsModel" begin
    systemData = MBD.TBPSystemData("Sun")
    dynamicsModel = MBD.TBPDynamicsModel(systemData)
    @test getStateSize(dynamicsModel, MBD.SIMPLE) == 6
    @test_throws ArgumentError appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1], MBD.SIMPLE)
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.SIMPLE) == [0.8234, 0, 0, 0, 0.1263, 0]
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL) == [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    @test evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.8234, 0, 0, 0, 0.1263, 0]) == [0, 0.1263, 0, -1.4749533162525872, 0, 0]
    @test getEquationsOfMotion(dynamicsModel, MBD.SIMPLE) == MBD.TBPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    @test isEpochIndependent(dynamicsModel)
    @test_throws ArgumentError getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [0, 0, 0, 0, 0, 0]
    @test_throws ArgumentError getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test size(getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])) == (6, 0)
    @test getPrimaryPosition(dynamicsModel) == [0.0, 0.0, 0.0]
    @test_throws ArgumentError getStateTransitionMatrix(dynamicsModel, [0.8324, 0, 0, 0, 0.1263, 0])
    @test getStateTransitionMatrix(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
    systemDataEarth = MBD.TBPSystemData("Earth")
    dynamicsModelEarth = MBD.TBPDynamicsModel(systemDataEarth)
    trajectory = getOsculatingOrbitalElements(dynamicsModel, [1.231789044172658E8, 0, 0, 0, 0.006519428065716839, 0])
    @test trajectory.E == -1077.3958252913806
end

@testset "TBPEquationsOfMotion" begin
    systemData = MBD.TBPSystemData("Sun")
    dynamicsModel = MBD.TBPDynamicsModel(systemData)
    EOMs_simple = MBD.TBPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    qdot_simple = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.SIMPLE))
    computeDerivatives!(qdot_simple, [0.8234, 0, 0, 0, 0.1263, 0], (EOMs_simple,), 0.0)
    @test qdot_simple == [0, 0.1263, 0, -1.4749533162525872, 0, 0]
    EOMs_full = MBD.TBPEquationsOfMotion(MBD.FULL, dynamicsModel)
    qdot_full = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.FULL))
    computeDerivatives!(qdot_full, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL), (EOMs_full,), 0.0)
    @test qdot_full == [0, 0.1263, 0, -1.4749533162525872, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 3.5825924611430353, 0, 0, 0, 0, 0, 0, -1.791296230571517, 0, 0, 0, 0, 0, 0, -1.791296230571517, 0, 0, 0]
    @test getLambertArc(dynamicsModel, [0.785798, 0.618484, -4.40907E-5], [-0.695636, -1.36224, -0.0114978], 5.391710131490492, "Long") == ([-61781.16548851388, 511386.43681523146, 7933.518512179655], [393204.0310118931, 137401.10858711743, -2458.8113571350514])
end

@testset "SPICE Functions" begin
    SPICE.furnsh("../src/spice/kernels/naif0012.tls", "../src/spice/kernels/de440.bsp", "../src/spice/kernels/mar097.bsp")
    EarthInitialState = getEphemerides("Nov 1 2026", [0.0], "Earth", "Sun", "ECLIPJ2000")
    @test EarthInitialState == [[1.1669364364177027E8, 9.184694649040842E7, -6547.627331614494, -18.896701708854742, 23.300520269960817, -0.0003837831634516675]]
    SPICE.kclear()
end

@testset "Utility Functions" begin
    @test_throws BoundsError checkIndices([1, 7], 6)
    @test_throws ArgumentError checkIndices([1, 1], 6)
    @test_throws ArgumentError maskData([true], [1.0 2.0 3.0; 4.0 5.0 6.0])
    mask1 = [true, false, true]
    data = [1.0 2.0 3.0; 4.0 5.0 6.0]
    @test maskData(mask1, data) == [1.0 3.0; 4.0 6.0]
    mask2 = [false, false, false]
    @test size(maskData(mask2, data)) == (1, 0)
    mask3 = [true, true, true]
    @test maskData(mask3, data) == data
end

@testset "Lyapunov Example" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    q0::Vector{Float64} = [0.8234, 0, 0, 0, 0.1263, 0]
    tSpan::Vector{Float64} = [0, 2.743]
    halfPeriod::Float64 = (tSpan[2]-tSpan[1])/2
    propagator = MBD.Propagator()
    arc::MBD.Arc = propagate(propagator, q0, tSpan, dynamicsModel)
    qf::Vector{Float64} = copy(getStateByIndex(arc, getStateCount(arc)))
    @test qf == [0.8322038365943337, -0.0027959784747644103, 0, 0.024609343452095516, 0.11660029449394844, 0]
    originNode = MBD.Node(tSpan[1], q0, dynamicsModel)
    originNode.state.name = "Initial State"
    terminalNode = MBD.Node(halfPeriod, qf, dynamicsModel)
    terminalNode.state.name = "Target State"
    segment = MBD.Segment(halfPeriod, originNode, terminalNode)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    multipleShooterProblem = MBD.MultipleShooterProblem()
    addSegment!(multipleShooterProblem, segment)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    x0Constraint = MBD.StateConstraint(originNode, [1], [q0[1]])
    qfConstraint = MBD.StateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(multipleShooterProblem, x0Constraint)
    addConstraint!(multipleShooterProblem, qfConstraint)
    addConstraint!(multipleShooterProblem, continuityConstraint)
    multipleShooter = MBD.MultipleShooter()
    solved::MBD.MultipleShooterProblem = solve!(multipleShooter, multipleShooterProblem)
    @test isapprox(solved.nodes[1].state.data-[0.8234, 0, 0, 0, 0.12623176201421427, 0], zeros(Float64, 6), atol = 1E-11)
end

@testset "Halo Multiple Shooter Example" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    nodeStates::Vector{Vector{Float64}} = [
        [1.0277, 0, 0.1857, 0, -0.1152, 0],
        [1.008, -0.0476, 0.1133, -0.0714, -0.0307, -0.2984],
        [1.008, 0.0476, 0.1133, 0.0714, -0.0307, 0.2985],
        [1.0277, 0, 0.1857, 0, -0.1152, 0]
    ]
    nodeTimes::Vector{Float64} = [0, 0.52854, 1.0571, 1.5856]
    nodes::Vector{MBD.Node} = Vector{MBD.Node}(undef, length(nodeTimes))
    [nodes[n] = MBD.Node(nodeTimes[n], nodeStates[n], dynamicsModel) for n in eachindex(nodeTimes)]
    multipleShooterProblem = MBD.MultipleShooterProblem()
    segments::Vector{MBD.Segment} = Vector{MBD.Segment}(undef, length(nodeTimes)-1)
    for s::Int64 in 1:length(nodeTimes)-1
        segments[s] = MBD.Segment(nodeTimes[s+1]-nodeTimes[s], nodes[s], nodes[s+1])
        addSegment!(multipleShooterProblem, segments[s])
    end
    @test length(multipleShooterProblem.segments) == 3
    map(s -> addConstraint!(multipleShooterProblem, MBD.ContinuityConstraint(s)), segments)
    addConstraint!(multipleShooterProblem, MBD.StateMatchConstraint(nodes[1].state, nodes[end].state, [1, 2, 3, 4, 5, 6]))
    addConstraint!(multipleShooterProblem, MBD.JacobiConstraint(nodes[1], 3.04))
    @test checkJacobian(multipleShooterProblem)
    multipleShooter = MBD.MultipleShooter()
    solved::MBD.MultipleShooterProblem = solve!(multipleShooter, multipleShooterProblem)
    @test isapprox(solved.nodes[1].state.data-solved.nodes[end].state.data, zeros(Float64, 6), atol = 1E-11)
end

@testset "Lyapunov Continuation Example" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    Moon = MBD.BodyData("Moon")
    L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
    (stateGuess::Vector{Float64}, timeGuess::Vector{Float64}) = getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])
    targetJC::Float64 = getJacobiConstant(dynamicsModel, stateGuess)
    targeter = LyapunovJCTargeter(dynamicsModel)
    solution1::MBD.MultipleShooterProblem = correct(targeter, stateGuess, timeGuess, targetJC)
    orbit1 = MBD.CR3BPPeriodicOrbit(solution1, targeter)
    solution2::MBD.MultipleShooterProblem = correct(targeter, stateGuess, timeGuess, 3.186)
    orbit2 = MBD.CR3BPPeriodicOrbit(solution2, targeter)
    continuationEngine = MBD.JacobiConstantContinuationEngine(-1E-3, -1E-2)
    continuationEngine.printProgress = false
    ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [NaN NaN; -2.5 0])
    addJumpCheck!(continuationEngine, ydot0JumpCheck)
    @test length(continuationEngine.jumpChecks) == 1
    numberSteps::Int64 = 100
    numberEndCheck = MBD.NumberStepsContinuationEndCheck(numberSteps)
    MoonEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [L1[1] getPrimaryPosition(dynamicsModel, 2)[1]-Moon.bodyRadius/systemData.charLength; NaN NaN])
    addEndCheck!(continuationEngine, numberEndCheck)
    addEndCheck!(continuationEngine, MoonEndCheck)
    @test length(continuationEngine.endChecks) == 2
    solutions::Vector{MBD.MultipleShooterProblem} = doContinuation!(continuationEngine, solution1, solution2)
    @test length(solutions) == numberSteps
end
