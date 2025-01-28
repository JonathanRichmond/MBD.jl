"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 9/1/22
U: 1/27/25
"""

using MBD, DifferentialEquations, LinearAlgebra, SPICE, StaticArrays
using Test

# include("ExampleBifurcations.jl")
include("ExampleLyapunovJCTargeter.jl")
include("ExampleLyapunovx0Targeter.jl")

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
    EarthData2 = MBD.BodyData("Earth")
    @test EarthData == EarthData2
    @test_throws ArgumentError MBD.CR3BPSystemData("Earth", "Mars")
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    @test systemData.primarySpiceIDs == SVector{2, Int16}([399, 301])
    @test systemData.primaryData[1] == EarthData
    systemData2 = MBD.CR3BPSystemData("Earth", "Moon")
    @test systemData == systemData2
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    @test dynamicsModel.systemData == systemData
    dynamicsModel2 = MBD.CR3BPDynamicsModel(systemData2)
    @test dynamicsModel == dynamicsModel2
    EOMs = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    @test EOMs.dynamicsModel == dynamicsModel
    EOMs2 = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel2)
    @test EOMs == EOMs2
    integratorFactory = MBD.IntegratorFactory(MBD.DP8)
    @test integratorFactory.integrator == DP8()
    integratorFactory2 = MBD.IntegratorFactory(MBD.DP8)
    @test integratorFactory == integratorFactory2
    propagator = MBD.Propagator(MBD.DP8, MBD.FULL)
    @test propagator.equationType == MBD.FULL
    @test propagator.integratorFactory == integratorFactory
    propagator2 = MBD.Propagator(MBD.DP8, MBD.FULL)
    @test propagator == propagator2
    arc = MBD.CR3BPArc(dynamicsModel)
    @test arc.dynamicsModel == dynamicsModel
    arc2 = MBD.CR3BPArc(dynamicsModel2)
    @test arc == arc2
    @test_throws ArgumentError MBD.Variable([1.0, 1.0], [true, false, true])
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test stateVariable.freeVariableMask == [true, false, false, false, true, false]
    stateVariable2 = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test stateVariable == stateVariable2
    @test_throws ArgumentError MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263], dynamicsModel)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    @test originNode.state.data == stateVariable.data
    originNode2 = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel2)
    @test originNode == originNode2
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    @test_throws ArgumentError MBD.CR3BPSegment(2.743, originNode, originNode)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    @test segment.terminalNode == terminalNode
    segment2 = MBD.CR3BPSegment(2.743, originNode2, terminalNode)
    @test segment == segment2
    problem = MBD.CR3BPMultipleShooterProblem()
    @test !problem.hasBeenBuilt
    problem2 = MBD.CR3BPMultipleShooterProblem()
    @test problem == problem2
    @test MBD.CR3BPContinuityConstraint <: MBD.AbstractConstraint
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    @test continuityConstraint.segment == segment
    continuityConstraint2 = MBD.CR3BPContinuityConstraint(segment2)
    @test continuityConstraint == continuityConstraint2
    @test MBD.CR3BPStateConstraint <: MBD.AbstractConstraint
    @test_throws ArgumentError MBD.CR3BPStateConstraint(originNode, [1, 5], [1.0])
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    @test_throws ArgumentError MBD.CR3BPStateConstraint(originNode, [2], [1.0])
    stateConstraint = MBD.CR3BPStateConstraint(originNode, [1, 5], [0.8, 0.1])
    @test stateConstraint.variable == originNode.state
    setFreeVariableMask!(originNode2.state, [true, false, false, false, true, false])
    stateConstraint2 = MBD.CR3BPStateConstraint(originNode2, [1, 5], [0.8, 0.1])
    @test stateConstraint == stateConstraint2
    @test MBD.StateMatchConstraint <: MBD.AbstractConstraint
    variable = MBD.Variable([1.0, 2.0], [true, false])
    setFreeVariableMask!(terminalNode.state, [true, false, true, true, true, true])
    @test_throws ArgumentError MBD.StateMatchConstraint(stateVariable, variable, [1, 2])
    @test_throws ArgumentError MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2])
    setFreeVariableMask!(terminalNode.state, [true, true, true, true, true, true])
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test stateMatchConstraint.variable1 == originNode.state
    stateMatchConstraint2 = MBD.StateMatchConstraint(originNode2.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test stateMatchConstraint == stateMatchConstraint2
    @test MBD.JacobiConstraint <: MBD.AbstractConstraint
    # systemDataTBP = MBD.TBPSystemData("Earth")
    # dynamicsModelTBP = MBD.TBPDynamicsModel(systemDataTBP)
    # nodeTBP = MBD.Node(0.0, originNode.state.data, dynamicsModelTBP)
    # @test_throws ArgumentError MBD.JacobiConstraint(nodeTBP, 3.0)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    @test JCConstraint.state == originNode.state
    JCConstraint2 = MBD.JacobiConstraint(originNode2, 3.0)
    @test JCConstraint == JCConstraint2
    convergenceCheck = MBD.ConstraintVectorL2NormConvergenceCheck()
    @test convergenceCheck.maxVectorNorm == 1E-10
    convergenceCheck2 = MBD.ConstraintVectorL2NormConvergenceCheck()
    @test convergenceCheck == convergenceCheck2
    @test MBD.MinimumNormUpdateGenerator <: MBD.AbstractUpdateGenerator
    @test MBD.LeastSquaresUpdateGenerator <: MBD.AbstractUpdateGenerator
    multipleShooter = MBD.CR3BPMultipleShooter()
    @test multipleShooter.maxIterations == Int16(25)
    multipleShooter2 = MBD.CR3BPMultipleShooter()
    @test multipleShooter == multipleShooter2
    continuationFamily = MBD.CR3BPContinuationFamily(dynamicsModel)
    @test continuationFamily.dynamicsModel == dynamicsModel
    continuationFamily2 = MBD.CR3BPContinuationFamily(dynamicsModel2)
    @test continuationFamily == continuationFamily2
    @test_throws ArgumentError MBD.AdaptiveStepSizeByElementGenerator("Initial State", -1, 1E-4, 1E-2)
    @test_throws ArgumentError MBD.AdaptiveStepSizeByElementGenerator("Initial State", 1, -1E-4, 1E-2)
    adaptiveStepSizeByElementGenerator = MBD.AdaptiveStepSizeByElementGenerator("Initial State", 1, -1E-4, -1E-2)
    @test adaptiveStepSizeByElementGenerator.elementIndex == Int16(1)
    adaptiveStepSizeByElementGenerator2 = MBD.AdaptiveStepSizeByElementGenerator("Initial State", 1, -1E-4, -1E-2)
    @test adaptiveStepSizeByElementGenerator == adaptiveStepSizeByElementGenerator2
    @test MBD.NumberStepsContinuationEndCheck <: MBD.AbstractContinuationEndCheck
    numberStepsContinuationEndCheck = MBD.NumberStepsContinuationEndCheck(100)
    @test numberStepsContinuationEndCheck.maxSteps == Int16(100)
    numberStepsContinuationEndCheck2 = MBD.NumberStepsContinuationEndCheck(100)
    @test numberStepsContinuationEndCheck == numberStepsContinuationEndCheck2
    @test MBD.BoundingBoxContinuationEndCheck <: MBD.AbstractContinuationEndCheck
    boundingBoxContinuationEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [0 0.99; NaN NaN])
    @test isequal(boundingBoxContinuationEndCheck.paramBounds, [0 0.99; NaN NaN])
    boundingBoxContinuationEndCheck2 = MBD.BoundingBoxContinuationEndCheck("Initial State", [0 0.99; NaN NaN])
    @test boundingBoxContinuationEndCheck == boundingBoxContinuationEndCheck2
    @test MBD.BoundingBoxJumpCheck <: MBD.AbstractContinuationJumpCheck
    boundingBoxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
    @test isequal(boundingBoxJumpCheck.paramBounds, [0 1.2; 0 1.0])
    boundingBoxJumpCheck2 = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
    @test boundingBoxJumpCheck == boundingBoxJumpCheck2
    # bifurcation = MBD.Bifurcation(orbitFamily, periodicOrbit, 1, MBD.TANGENT, 1)
    # @test bifurcation.orbit == 1
    # periodicOrbit = MBD.CR3BPPeriodicOrbit(problem, 0.0, zeros(Float64, (6,6))
    # @test periodicOrbit.dynamicsModel == dynamicsModel
    # periodicOrbit2 = MBD.CR3BPPeriodicOrbit(problem2, 0.0, zeros(Float64, (6,6))
    # @test periodicOrbit == periodicOrbit2
    # @test MBD.CR3BPOrbitFamily <: MBD.AbstractStructureFamily
    # orbitFamily = MBD.CR3BPOrbitFamily([periodicOrbit])
    # @test length(orbitFamily.familyMembers) == 1
    # @test MBD.CR3BPManifoldArc <: MBD.AbstractTrajectoryStructure
    # manifoldArc = MBD.CR3BPManifoldArc([0.8234+0im, 0, 0, 0, 0.1263+0im, 0], periodicOrbit, 0.0)
    # @test manifoldArc.TOF == 0
    # @test MBD.TBPSystemData <: MBD.AbstractSystemData
    # @test systemDataTBP.charTime == 2.898139774118648e9
    # @test MBD.TBPDynamicsModel <: MBD.AbstractDynamicsModel
    # @test dynamicsModelTBP.systemData == systemDataTBP
    # @test MBD.TBPEquationsOfMotion <: MBD.AbstractEquationsOfMotion
    # EOMsTBP = MBD.TBPEquationsOfMotion(MBD.SIMPLE, dynamicsModelTBP)
    # @test EOMsTBP.mu == systemDataTBP.gravParam
    # @test MBD.TBPTrajectory <: MBD.AbstractTrajectoryStructure
    # trajectory = MBD.TBPTrajectory(originNode.state.data, dynamicsModelTBP)
    # @test trajectory.E == 0
end

@testset "Copy" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    @test MBD.shallowClone(originNode) == originNode
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    @test MBD.shallowClone(segment) == segment
    problem = MBD.CR3BPMultipleShooterProblem()
    @test MBD.shallowClone(problem) == problem
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    @test MBD.shallowClone(continuityConstraint) == continuityConstraint
    stateConstraint = MBD.CR3BPStateConstraint(originNode, [1, 5], [0.8, 0.1])
    @test MBD.shallowClone(stateConstraint) == stateConstraint
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test MBD.shallowClone(stateMatchConstraint) == stateMatchConstraint
    jacobiConstraint = MBD.JacobiConstraint(originNode, 3.0)
    @test MBD.shallowClone(jacobiConstraint) == jacobiConstraint
#     systemDataTBP = MBD.TBPSystemData("Earth")
#     dynamicsModelTBP = MBD.TBPDynamicsModel(systemDataTBP)
#     trajectory = MBD.TBPTrajectory(originNode.state.data, dynamicsModelTBP)
#     @test shallowClone(trajectory) == trajectory
end

@testset "Deep Copy" begin
    variable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test MBD.deepClone(variable) == variable
    problem = MBD.CR3BPMultipleShooterProblem()
    @test MBD.deepClone(problem) == problem
end

@testset "BodyName" begin
    bodyName = MBD.BodyName("Earth")
    @test getIDCode(bodyName) == Int16(399)
end

@testset "CR3BPSystemData" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    @test isApproxSigFigs([getCharLength(systemData)], [384747.9920112924], 13)
    @test isApproxSigFigs([getCharMass(systemData)], [6.046042990276359E24], 13)
    @test isApproxSigFigs([getCharTime(systemData)], [375699.8590849912], 13)
    @test isApproxSigFigs([getMassRatio(systemData)], [0.012150584269940356], 13)
    @test getNumPrimaries(systemData) == Int16(2)
end

@testset "CR3BPDynamicsModel" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    @test_throws ArgumentError appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1.0], MBD.SIMPLE)
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.SIMPLE) == [0.8234, 0, 0, 0, 0.1263, 0]
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL) == [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.ARCLENGTH) == [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0]
    @test isApproxSigFigs(evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.8234, 0, 0, 0, 0.1263, 0]), [0, 0.1263, 0, 0.1103323864940, 0, 0], 13)
    @test isApproxSigFigs(evaluateEquations(dynamicsModel, MBD.FULL, 0.0, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL)), [0, 0.1263, 0, 0.1103323864940, 0, 0, 0, 0, 0, 9.851145859594, 0, 0, 0, 0, 0, 0, -3.425572929797, 0, 0, 0, 0, 0, 0, -4.425572929797, 1.0, 0, 0, 0, -2.0, 0, 0, 1.0, 0, 2.0, 0, 0, 0, 0, 1.0, 0, 0, 0], 13)
    @test isApproxSigFigs(evaluateEquations(dynamicsModel, MBD.ARCLENGTH, 0.0, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.ARCLENGTH)), [0, 0.1263, 0, 0.1103323864940, 0, 0, 0, 0, 0, 9.851145859594, 0, 0, 0, 0, 0, 0, -3.425572929797, 0, 0, 0, 0, 0, 0, -4.425572929797, 1.0, 0, 0, 0, -2.0, 0, 0, 1.0, 0, 2.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0.1263], 13)
    MoonData = MBD.BodyData("Moon")
    @test isApproxSigFigs(get2BApproximation(dynamicsModel, MoonData, 2, 0.01), [0.9778494157301, 0, 0, 0, 1.112296886957, 0], 13)
    @test isApproxSigFigs([getCharLength(dynamicsModel)], [384747.9920112924], 13)
    @test isApproxSigFigs([getCharTime(dynamicsModel)], [375699.8590849912], 13)
    @test getEquationsOfMotion(dynamicsModel, MBD.SIMPLE) == MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    @test_throws ArgumentError getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [0, 0, 0, 0, 0, 0]
    @test_throws ArgumentError getEquilibriumPoint(dynamicsModel, 6)
    @test getEquilibriumPoint(dynamicsModel, 1) == [0.8369151323643023, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, 2) == [1.1556821602923404, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, 3) == [-1.005062645252109, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, 4) == [0.48784941573005963, 0.8660254037844386, 0]
    @test isApproxSigFigs([getExcursion(dynamicsModel, 2, [0.8234, 0, 0, 0, 0.1263, 0])], [63271.58248957], 13)
    @test isApproxSigFigs([getJacobiConstant(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])], [3.174356023206], 13)
    L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
    L5::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 5)
    @test_throws ArgumentError getLinearVariation(dynamicsModel, 6, L1, [0.005, 0, 0])
    @test_throws ArgumentError getLinearVariation(dynamicsModel, 5, L5, [0.005, 0, 0], "Medium")
    @test isApproxSigFigs(getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])[1], [0.8419151323643, 0, 0, 0, -0.04186136597443, 0], 13)
    @test isApproxSigFigs(getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])[2], [0, 2.691579560798], 13)
    @test isApproxSigFigs(getLinearVariation(dynamicsModel, 5, L5, [0.005, 0, 0], "Long")[1], [0.4928494157301, -0.8660254037844, 0, -0.003168674904327, -0.002097320259364, 0], 13)
    @test isApproxSigFigs(getLinearVariation(dynamicsModel, 5, L5, [0.005, 0, 0], "Long")[2], [0, 21.06979705456], 13)
    @test isApproxSigFigs([getMassRatio(dynamicsModel)], [0.012150584269940356], 13)
    @test_throws ArgumentError getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test size(getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])) == (6, 0)
    @test_throws ArgumentError getPrimaryPosition(dynamicsModel, 3)
    @test isApproxSigFigs(getPrimaryPosition(dynamicsModel, 1), [-0.01215058426994, 0, 0], 13)
    @test isApproxSigFigs(getPseudopotentialJacobian(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0]), [9.851145859594, -3.425572929797, -4.425572929797, 0, 0, 0], 13)
    @test getStateSize(dynamicsModel, MBD.ARCLENGTH) == Int16(43)
    @test_throws ArgumentError getStateTransitionMatrix(dynamicsModel, [0.8324, 0, 0, 0, 0.1263, 0])
    @test getStateTransitionMatrix(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0]) == [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
    @test_throws ArgumentError getTidalAcceleration(dynamicsModel, 3, [0.8234, 0, 0])
    @test isApproxSigFigs(getTidalAcceleration(dynamicsModel, 1, [0.8234, 0, 0]), [-0.5915635514335, 0, 0], 13)
    @test isApproxSigFigs(getTidalAcceleration(dynamicsModel, 2, [0.8234, 0, 0]), [1.272695937927, 0, 0], 13)
    @test isEpochIndependent(dynamicsModel)
    @test_throws ArgumentError primaryInertial2Rotating(dynamicsModel, 2, [[1.0, 0, 0, 0, 1.0, 0], [1.1, 0, 0, 0, 1.1, 0]], [0.0])
    @test_throws ArgumentError primaryInertial2Rotating(dynamicsModel, 3, [[1.0, 0, 0, 0, 1.0, 0]], [0.0])
    @test isApproxSigFigs(primaryInertial2Rotating(dynamicsModel, 1, [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[1], [0.8112494157301, 0, 0, 0, -0.6971, 0], 13)
    @test isApproxSigFigs(primaryInertial2Rotating(dynamicsModel, 1, [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[2], [-0.7802015843359, -0.3204195756524, 0, -0.2978446415273, 0.6510398005646, 0], 13)
    SPICE.furnsh("../src/spice/kernels/naif0012.tls", "../src/spice/kernels/de440.bsp", "../src/spice/kernels/mar097.bsp")
    @test isApproxSigFigs(rotating2PrimaryEclipJ2000(dynamicsModel, "Jan 1 2026", [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[1], [0.3093978367126, 0.7726559358880, 0.07362446772986, -0.8930326868322, 0.3572626063404, 0.003549199984676], 13)
    @test isApproxSigFigs(rotating2PrimaryEclipJ2000(dynamicsModel, "Jan 1 2026", [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[2], [-0.5944039158312, -0.5959053293302, -0.06725346738315, 0.6649246782108, -0.6932473627690, -0.03831818053376], 13)
    @test isApproxSigFigs(rotating2PrimaryInertial(dynamicsModel, 1, [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[1], [0.8355505842699, 0, 0, 0, 0.9618505842699, 0], 13)
    @test isApproxSigFigs(rotating2PrimaryInertial(dynamicsModel, 1, [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[2], [-0.7770787174211, 0.3302890741241, 0, -0.3982243502394, -0.8749870761992, 0], 13)
    @test isApproxSigFigs(rotating2SunEclipJ2000(dynamicsModel, "Jan 1 2026", [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[1], [0.3338038885097, 0.7624427656725, 0.07348994464091, -0.8813086468342, 0.3852438101716, 0.006230746360073], 13)
    @test isApproxSigFigs(rotating2SunEclipJ2000(dynamicsModel, "Jan 1 2026", [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[2], [-0.6130761226760, -0.5767984348495, -0.06620754424084, 0.6426274696210, -0.7138340110843, -0.04069344599638], 13)
    @test isApproxSigFigs(secondaryEclipJ20002Rotating(dynamicsModel, "Jan 1 2026", [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[1], [1.292747986260, -0.7644878803039, -0.02413923274187, -0.6476948978804, -0.2640805003305, -0.01050524106962], 13)
    @test isApproxSigFigs(secondaryEclipJ20002Rotating(dynamicsModel, "Jan 1 2026", [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[2], [0.4059409161323, 0.5944489928956, -0.02416476991132, 0.4946209535941, 0.5115726248560, -0.01041991035573], 13)
    SPICE.kclear()
end

@testset "CR3BPEquationsOfMotion" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    EOMs_simple = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    qdot_simple = Vector{Float64}(undef, getStateSize(EOMs_simple))
    computeDerivatives!(qdot_simple, [0.8234, 0, 0, 0, 0.1263, 0], (EOMs_simple,), 0.0)
    @test isApproxSigFigs(qdot_simple, [0, 0.1263, 0, 0.1103323864940, 0, 0], 13)
    EOMs_full = MBD.CR3BPEquationsOfMotion(MBD.FULL, dynamicsModel)
    qdot_full = Vector{Float64}(undef, getStateSize(EOMs_full))
    computeDerivatives!(qdot_full, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL), (EOMs_full,), 0.0)
    @test isApproxSigFigs(qdot_full, [0, 0.1263, 0, 0.1103323864940, 0, 0, 0, 0, 0, 9.851145859594, 0, 0, 0, 0, 0, 0, -3.425572929797, 0, 0, 0, 0, 0, 0, -4.425572929797, 1, 0, 0, 0, -2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0], 13)
    EOMs_arclength = MBD.CR3BPEquationsOfMotion(MBD.ARCLENGTH, dynamicsModel)
    qdot_arclength = Vector{Float64}(undef, getStateSize(EOMs_arclength))
    computeDerivatives!(qdot_arclength, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.ARCLENGTH), (EOMs_arclength,), 0.0)
    @test isApproxSigFigs(qdot_arclength, [0, 0.1263, 0, 0.1103323864940, 0, 0, 0, 0, 0, 9.851145859594, 0, 0, 0, 0, 0, 0, -3.425572929797, 0, 0, 0, 0, 0, 0, -4.425572929797, 1, 0, 0, 0, -2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0.1263], 13)
    EOMs_momentum = MBD.CR3BPEquationsOfMotion(MBD.MOMENTUM, dynamicsModel)
    qdot_momentum = Vector{Float64}(undef, getStateSize(EOMs_momentum))
    computeDerivatives!(qdot_momentum, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.MOMENTUM), (EOMs_momentum,), 0.0)
    @test isApproxSigFigs(qdot_momentum, [0, 0.1263, 0, 0.1103323864940, 0, 0, 0, 0, 0, 9.851145859594, 0, 0, 0, 0, 0, 0, -3.425572929797, 0, 0, 0, 0, 0, 0, -4.425572929797, 1, 0, 0, 0, -2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0], 13)
    @test isApproxSigFigs([getMassRatio(EOMs_simple)], [0.012150584269940356], 13)
end

@testset "Propagator" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.CR3BPArc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test isApproxSigFigs(getStateByIndex(arc, -1), [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], 10)
    @test isApproxSigFigs([getTimeByIndex(arc, -2)], [2.685938962364293], 10)
    distanceEvent = DifferentialEquations.ContinuousCallback(p2DistanceCondition, terminateAffect!)
    eventArc::MBD.CR3BPArc = propagateWithEvent(propagator, distanceEvent, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel, [getPrimaryPosition(dynamicsModel, 2)[1]-0.83])
    @test isApproxSigFigs(getStateByIndex(eventArc, -1), [0.8399883237735, 0.05525880501532, 0, 0.03871870586244, 0.01970130539743, 0], 10)
end

@testset "CR3BPArc" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.CR3BPArc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test isApproxSigFigs([getMassRatio(arc)], [0.012150584269940356], 13)
    stateCount::Int64 = getStateCount(arc)
    @test_throws BoundsError getStateByIndex(arc, stateCount+1)
    @test isApproxSigFigs(getStateByIndex(arc, 20), [0.8539411650664558, -0.029778764729793614, 0, -0.009955443086921378, -0.11231429456648163, 0], 10)
    @test isApproxSigFigs(getStateByIndex(arc, -1), [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], 10)
    @test_throws BoundsError deleteStateAndTime!(arc, stateCount+1)
    deleteStateAndTime!(arc, stateCount)
    @test_throws BoundsError getTimeByIndex(arc, stateCount)
    @test isApproxSigFigs([getTimeByIndex(arc, 20)], [1.6057680454134118], 10)
    @test isApproxSigFigs([getTimeByIndex(arc, -1)], [2.685938962364293], 10)
end

@testset "Variable" begin
    variable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test getData(variable) == [0.8234, 0, 0, 0, 0.1263, 0]
    @test getFreeVariableData(variable) == [0.8234, 0.1263]
    @test getFreeVariableMask(variable) == [true, false, false, false, true, false]
    @test getNumFreeVariables(variable) == 2
    @test_throws ArgumentError setFreeVariableData!(variable, [4.0, 7.0, 8.5])
    setFreeVariableData!(variable, [4.0, 7.0])
    @test variable.data == [4.0, 0, 0, 0, 7.0, 0]
    @test_throws ArgumentError setFreeVariableMask!(variable, [true, false])
    setFreeVariableMask!(variable, [true, false, false, false, false, false])
    @test variable.freeVariableMask == [true, false, false, false, false, false]
end

@testset "CR3BPNode" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    node = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, true, true, true, true, true])
    stateVariable.name = "Node State"
    epochVariable = MBD.Variable([0.0], [false])
    epochVariable.name = "Node Epoch"
    @test getVariables(node) == [stateVariable, epochVariable]
end

@testset "CR3BPSegment" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    @test isApproxSigFigs(getFinalState!(segment), [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], 10)
    @test isApproxSigFigs(getFinalStateRate!(segment), [0.024609343451907937, 0.11660029449401914, 0, 0.18113477226039054, -0.03842090358792116, 0], 10)
    @test getPartials_FinalStateWRTEpoch!(segment) == zeros(Float64, (6,1))
    @test isApproxSigFigs(getPartials_FinalStateWRTInitialState!(segment)[1,:], [1363.7623463360135, -350.51893299267806, 0, 401.47006684119657, 130.91187990878421, 0], 10)
    TOFVariable = MBD.Variable([2.743], [true])
    TOFVariable.name = "Segment Time-of-Flight"
    @test getVariables(segment) == [TOFVariable]
    lazyPropagate!(segment, MBD.FULL)
    @test isApproxSigFigs(getStateByIndex(segment.propArc, -1), [0.8322038366125994, -0.0027959784802348967, 0, 0.024609343504035142, 0.11660029447417915, 0, 1363.7623463360135, -421.0811592085723, 0, 3877.9818867295958, -1465.21468708984, 0, -350.51893299267806, 109.07625516691418, 0, -995.8827212582189, 376.8074206183672, 0, 0, 0, 0.9853225289766858, 0, 0, -0.09704680343478124, 401.47006684119657, -123.80526849072281, 0, 1141.6513905844304, -431.6882797633952, 0, 130.91187990878421, -40.55137319725648, 0, 372.4229726003141, -139.7310684169993, 0, 0, 0, -0.07418084339358999, 0, 0, 1.0222023592351468], 10)
    resetPropagatedArc!(segment)
    @test getStateCount(segment.propArc) == 1
end

@testset "CR3BPMultipleShooterProblem" begin
    problem = MBD.CR3BPMultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    @test length(problem.segments) == 1
    @test_throws ErrorException buildAdjacencyMatrix(problem)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    addConstraint!(problem, continuityConstraint)
    @test length(problem.constraintIndexMap) == 1
    badProblem = MBD.CR3BPMultipleShooterProblem()
    stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    addVariable!(badProblem, stateVariable)
    @test length(badProblem.freeVariableIndexMap) == 1
    buildProblem!(problem)
    @test problem.hasBeenBuilt
    @test checkJacobian(problem)
    @test length(getConstraints(problem)) == 1
    @test isApproxSigFigs(getConstraintVector!(problem), [1.4937372928680581e-5, -4.625045526685141e-6, -1.0e-8, 4.249830200013863e-5, -1.6059545274130227e-5, -1.0e-8], 10)
    @test length(getFreeVariableIndexMap!(problem)) == 5
    @test isApproxSigFigs([LinearAlgebra.norm(getFreeVariableVector!(problem))], [2.987433787528257], 10)
    @test size(getJacobian!(problem)) == (6, 9)
    @test isApproxSigFigs([LinearAlgebra.norm(getJacobian!(problem)[1,:])], [1370.1021037432147], 10)
    @test getNumConstraints(problem) == 6
    @test getNumFreeVariables!(problem) == 9
    removeConstraint!(problem, collect(keys(problem.constraintIndexMap))[1])
    @test length(problem.constraintIndexMap) == 0
    newFreeVariableVector::Vector{Float64} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    setFreeVariableVector!(problem, newFreeVariableVector)
    @test getFreeVariableVector!(problem) == newFreeVariableVector
end

@testset "CR3BPContinuityConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    problem = MBD.CR3BPMultipleShooterProblem()
    @test isApproxSigFigs(evaluateConstraint(continuityConstraint, problem.freeVariableIndexMap, problem.freeVariableVector), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 10)
    @test getNumConstraintRows(continuityConstraint) == 6
    partials::Dict{MBD.Variable, Matrix{Float64}} = getPartials_ConstraintWRTVariables(continuityConstraint, problem.freeVariableIndexMap, problem.freeVariableVector)
    @test length(keys(partials)) == 3
    @test isApproxSigFigs(partials[segment.TOF][:,1], [0.024609343504035142, 0.11660029447417915, 0, 0.18113477241637266, -0.03842090366810482, 0], 10)
end

@testset "CR3BPStateConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    stateConstraint = MBD.CR3BPStateConstraint(originNode, [1], [originNode.state.data[1]])
    problem = MBD.CR3BPMultipleShooterProblem()
    @test evaluateConstraint(stateConstraint, problem.freeVariableIndexMap, problem.freeVariableVector) == [0]
    @test getNumConstraintRows(stateConstraint) == 1
    partials::Dict{MBD.Variable, Matrix{Float64}} = getPartials_ConstraintWRTVariables(stateConstraint, problem.freeVariableIndexMap, problem.freeVariableVector)
    @test length(keys(partials)) == 1
    @test isApproxSigFigs(partials[originNode.state][1,:], [1.0, 0, 0, 0, 0, 0], 10)
end

@testset "StateMatchConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    problem = MBD.CR3BPMultipleShooterProblem()
    @test isApproxSigFigs(evaluateConstraint(stateMatchConstraint, problem.freeVariableIndexMap, problem.freeVariableVector), [-0.008803836594267866, 0.0027959784747437957, 0, -0.024609343451907937, 0.009699705505980857, 0], 10)
    @test getNumConstraintRows(stateMatchConstraint) == 6
    partials::Dict{MBD.Variable, Matrix{Float64}} = getPartials_ConstraintWRTVariables(stateMatchConstraint, problem.freeVariableIndexMap, problem.freeVariableVector)
    @test length(keys(partials)) == 2
    @test isApproxSigFigs(partials[originNode.state][:,1], [1.0, 0, 0, 0, 0, 0], 10)
end

@testset "JacobiConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.1743)
    problem = MBD.CR3BPMultipleShooterProblem()
    @test isApproxSigFigs(evaluateConstraint(JCConstraint, problem.freeVariableIndexMap, problem.freeVariableVector), [5.6023205926347686e-5], 10)
    @test getNumConstraintRows(JCConstraint) == 1
    partials::Dict{MBD.Variable, Matrix{Float64}} = getPartials_ConstraintWRTVariables(JCConstraint, problem.freeVariableIndexMap, problem.freeVariableVector)
    @test length(keys(partials)) == 1
    @test isApproxSigFigs(partials[JCConstraint.state][1,:], [-0.2845352270120183, 0, 0, 0, -0.2526, 0], 10)
end

@testset "ConstraintVectorL2NormConvergenceCheck" begin
    convergenceCheck = MBD.ConstraintVectorL2NormConvergenceCheck()
    problem = MBD.CR3BPMultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    addConstraint!(problem, continuityConstraint)
    @test isConverged(convergenceCheck, problem)
end

@testset "MinimumNormUpdateGenerator" begin
    minimumNormUpdate = MBD.MinimumNormUpdateGenerator()
    problem = MBD.CR3BPMultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    stateConstraint = MBD.CR3BPStateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(problem, continuityConstraint)
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, stateConstraint)
    @test canGenerateUpdate(minimumNormUpdate, problem)
    @test isApproxSigFigs([LinearAlgebra.norm(getFullUpdate(minimumNormUpdate, problem))], [1.282003911813496], 8)
    initialStateConstraint = MBD.CR3BPStateConstraint(originNode, [1], [0.82])
    addConstraint!(problem, initialStateConstraint)
    @test !canGenerateUpdate(minimumNormUpdate, problem)
    @test_throws ErrorException getFullUpdate(minimumNormUpdate, problem)
end

@testset "LeastSquaresUpdateGenerator" begin
    leastSquaresUpdate = MBD.LeastSquaresUpdateGenerator()
    problem = MBD.CR3BPMultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    addSegment!(problem, segment)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    addConstraint!(problem, continuityConstraint)
    @test !canGenerateUpdate(leastSquaresUpdate, problem)
    @test_throws ErrorException getFullUpdate(leastSquaresUpdate, problem)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    stateConstraint = MBD.CR3BPStateConstraint(terminalNode, [2, 4, 6], [0.0, 0.0, 0.0])
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, stateConstraint)
    @test canGenerateUpdate(leastSquaresUpdate, problem)
    @test isApproxSigFigs([LinearAlgebra.norm(getFullUpdate(leastSquaresUpdate, problem))], [1.2820039115787545], 8)
end

@testset "CR3BPMultipleShooter" begin
    multipleShooter = MBD.CR3BPMultipleShooter()
    badProblem = MBD.CR3BPMultipleShooterProblem()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    addSegment!(badProblem, segment)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.0)
    stateConstraint = MBD.CR3BPStateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(badProblem, continuityConstraint)
    addConstraint!(badProblem, JCConstraint)
    addConstraint!(badProblem, stateConstraint)
    @test_throws ErrorException MBD.solve!(multipleShooter, badProblem)
    problem = MBD.CR3BPMultipleShooterProblem()
    addSegment!(problem, segment)
    newJCConstraint = MBD.JacobiConstraint(originNode, 3.1743)
    addConstraint!(problem, continuityConstraint)
    addConstraint!(problem, newJCConstraint)
    addConstraint!(problem, stateConstraint)
    @test isApproxSigFigs([LinearAlgebra.norm(solveUpdateEquation(multipleShooter, problem))], [0.02797361726725746], 8)
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(multipleShooter, problem)
    @test isConverged(multipleShooter.convergenceCheck, solution)
end

@testset "CR3BPContinuationFamily" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    continuationFamily = MBD.CR3BPContinuationFamily(dynamicsModel)
    @test getNumMembers(continuationFamily) == 0
end

# @testset "CR3BPContinuationData" begin
# end

# @testset "AdaptiveStepSizeByElementGenerator" begin
#     stepSizeGenerator = MBD.AdaptiveStepSizeByElementGenerator("Initial State", 1, -1E-4, -1E-2)
#     continuationData = MBD.ContinuationData()
#     continuationData.currentStepSize = stepSizeGenerator.initialStepSize
#     updateStepSize!(stepSizeGenerator, continuationData)
#     @test continuationData.currentStepSize == -0.0002
#     continuationData.numIterations = 15
#     updateStepSize!(stepSizeGenerator, continuationData)
#     @test continuationData.currentStepSize == -0.0001
#     continuationData.numIterations = 1
#     continuationData.currentStepSize = -0.05
#     updateStepSize!(stepSizeGenerator, continuationData)
#     @test continuationData.currentStepSize == -0.01
#     continuationData.currentStepSize = stepSizeGenerator.initialStepSize
#     continuationData.converging = false
#     updateStepSize!(stepSizeGenerator, continuationData)
#     @test continuationData.currentStepSize == -0.00005
#     continuationData.currentStepSize = -1E-10
#     updateStepSize!(stepSizeGenerator, continuationData)
#     @test continuationData.forceEndContinuation
# end

# @testset "NumberStepsContinuationEndCheck" begin
#     stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
#     continuationData = MBD.ContinuationData()
#     @test !isContinuationDone(stepsEndCheck, continuationData)
#     continuationData.stepCount = 100
#     oldstd = stdout
#     redirect_stdout(devnull)
#     @test isContinuationDone(stepsEndCheck, continuationData)
#     redirect_stdout(oldstd)
# end

# @testset "BoundingBoxContinuationEndCheck" begin
#     boxEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [0 0.99; NaN NaN])
#     stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
#     stateVariable.name = "Initial State"
#     checkBounds(boxEndCheck, stateVariable)
#     setFreeVariableMask!(stateVariable, [true, false, false, false, false, false])
#     @test_throws ArgumentError checkBounds(boxEndCheck, stateVariable)
#     boxEndCheck.paramBounds = [0 0.5 0.99]
#     @test_throws ArgumentError checkBounds(boxEndCheck, stateVariable)
#     boxEndCheck.paramBounds = [0.5 0]
#     @test_throws ArgumentError checkBounds(boxEndCheck, stateVariable)
#     #isContinuationDone
# end

# @testset "BoundingBoxJumpCheck" begin
#     boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 0.99; NaN NaN])
#     problem = MBD.MultipleShooterProblem()
#     systemData = MBD.CR3BPSystemData("Earth", "Moon")
#     dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
#     originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
#     setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
#     terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
#     segment = MBD.Segment(2.743, originNode, terminalNode)
#     addSegment!(problem, segment)
#     addBounds!(boxJumpCheck, problem, problem.segments[1].originNode.state, [0 0.9; NaN NaN])
#     @test length(boxJumpCheck.variableBounds) == 1
#     stateVariable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
#     @test_throws MethodError addBounds!(boxJumpCheck, problem, stateVariable, [0 0.89; NaN NaN])
#     setFreeVariableMask!(stateVariable, [true, false, false, false, false, false])
#     @test_throws ArgumentError checkBounds(boxJumpCheck, stateVariable, [0 0.9; NaN NaN])
#     @test_throws ArgumentError checkBounds(boxJumpCheck, stateVariable, [0 0.5 0.9])
#     @test_throws ArgumentError checkBounds(boxJumpCheck, stateVariable, [0.9 0])
#     removeBounds!(boxJumpCheck, problem, problem.segments[1].originNode.state)
#     @test length(boxJumpCheck.variableBounds) == 0
#     @test_throws MethodError removeBounds!(boxJumpCheck, problem, stateVariable)
#     #isFamilyMember
# end

# @testset "NaturalParameterContinuationEngine" begin
#     naturalParamContinuation = MBD.NaturalParameterContinuationEngine("Initial State", 1, -1E-4, -1E-2)
#     stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
#     addEndCheck!(naturalParamContinuation, stepsEndCheck)
#     @test length(naturalParamContinuation.endChecks) == 1
#     boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
#     addJumpCheck!(naturalParamContinuation, boxJumpCheck)
#     @test length(naturalParamContinuation.jumpChecks) == 1
#     naturalParamContinuation.dataInProgress.converging = false
#     resetEngine!(naturalParamContinuation)
#     @test naturalParamContinuation.dataInProgress.converging
#     #computeFullStep
#     #constrainNextGuess!
#     #convergeInitialSolution
#     #doContinuation!
#     #endContinuation
#     #tryConverging!
# end

# @testset "JacobiConstantContinuationEngine" begin
#     JCContinuation = MBD.JacobiConstantContinuationEngine(-1E-4, -1E-1)
#     stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
#     addEndCheck!(JCContinuation, stepsEndCheck)
#     @test length(JCContinuation.endChecks) == 1
#     boxJumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0 1.2; 0 1.0])
#     addJumpCheck!(JCContinuation, boxJumpCheck)
#     @test length(JCContinuation.jumpChecks) == 1
#     JCContinuation.dataInProgress.converging = false
#     resetEngine!(JCContinuation)
#     @test JCContinuation.dataInProgress.converging
#     #computeFullStep
#     #constrainNextGuess!
#     #convergeInitialSolution
#     #doContinuation!
#     #endContinuation
#     #tryConverging!
# end

# @testset "Bifurcation" begin    
# end

@testset "CR3BPPeriodicOrbit and CR3BPManifold and CR3BPManifoldArc" begin
    shooter = MBD.CR3BPMultipleShooter()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.CR3BPNode(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false])
    terminalNode = MBD.CR3BPNode(2.743, [0.8322038365942679, -0.0027959784747437957, 0, 0.024609343451907937, 0.11660029449401914, 0], dynamicsModel)
    segment = MBD.CR3BPSegment(2.743, originNode, terminalNode)
    problem = MBD.CR3BPMultipleShooterProblem()
    addSegment!(problem, segment)
    continuityConstraint = MBD.CR3BPContinuityConstraint(segment)
    JCConstraint = MBD.JacobiConstraint(originNode, 3.1743)
    stateConstraint = MBD.CR3BPStateConstraint(terminalNode, [2, 4], [0.0, 0.0])
    addConstraint!(problem, continuityConstraint)
    addConstraint!(problem, JCConstraint)
    addConstraint!(problem, stateConstraint)
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)
    period::Float64 = 2*solution.segments[1].TOF.data[1]
    propagator = MBD.Propagator()
    propagator.equationType = MBD.STM
    initialState::Vector{Float64} = appendExtraInitialConditions(dynamicsModel, solution.nodes[1].state.data, MBD.STM)
    arc::MBD.CR3BPArc = propagate(propagator, initialState, [0, period/2], dynamicsModel)
    stateSize::Int64 = getStateSize(dynamicsModel, MBD.STM)
    halfPeriodState::StaticArrays.SVector{stateSize, Float64} = StaticArrays.SVector{stateSize, Float64}(getStateByIndex(arc, -1))
    halfPeriodSTM::StaticArrays.SMatrix{6, 6, Float64} = [halfPeriodState[7:12] halfPeriodState[13:18] halfPeriodState[19:24] halfPeriodState[25:30] halfPeriodState[31:36] halfPeriodState[37:42]]
    G::StaticArrays.SMatrix{6, 6, Float64} = SMatrix{6, 6}(1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, -1.0)
    Omega::StaticArrays.SMatrix{3, 3, Float64} = SMatrix{3, 3, }(0, -1.0, 0, 1.0, 0, 0, 0, 0, 0)
    A::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([zeros(Float64, (3,3)) [-1.0 0 0; 0 -1.0 0; 0 0 -1.0]; [1.0 0 0; 0 1.0 0; 0 0 1.0] -2*Omega])
    B::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([-2*Omega [1.0 0 0; 0 1.0 0; 0 0 1.0]; [-1.0 0 0; 0 -1.0 0; 0 0 -1.0] zeros(Float64, (3,3))])
    monodromy::Matrix{Float64} = G*A*(halfPeriodSTM')*B*G*halfPeriodSTM
    periodicOrbit = MBD.CR3BPPeriodicOrbit(solution.nodes[1].dynamicsModel, solution.nodes[1].state.data[1:6], period, monodromy)
    periodicOrbit2 = MBD.CR3BPPeriodicOrbit(solution.nodes[1].dynamicsModel, solution.nodes[1].state.data[1:6], period, monodromy)
    @test periodicOrbit == periodicOrbit2
    @test isApproxSigFigs(getBrouckeStability(periodicOrbit), [-5.569781483265836e6, 1.11422045234375E7], 8)
    @test isapprox(real.(getEigenData(periodicOrbit)[1]), [1.7994172955675235E-7, 0.9784501434972475, 0.9999942784866074, 1.0000057216849134, 1.0220244809057595, 5.569779482791029E6])
    @test isApproxSigFigs([getJacobiConstant(periodicOrbit)], [3.1743], 8)
    d::Float64 = 25/getCharLength(systemData)
    manifoldArc::MBD.CR3BPManifoldArc = getManifoldArcByTime(periodicOrbit, "Unstable", "Positive", d, 0.5)
    @test isapprox(real.(manifoldArc.initialCondition), [0.8234305489801728, -2.015402439111196E-5, 0, 0.00016950088194608143, 0.12648714319946022, 0])
    UPManifold::MBD.CR3BPManifold = getManifoldByArclength(periodicOrbit, "Unstable", "Positive", d, 10)
    @test isapprox(real.(UPManifold.initialConditions[2]), [0.8541164517830968, 0.024428003606702244, 0, 0.008747748408137988, -0.1195384951544451, 0])
    UNManifold::MBD.CR3BPManifold = getManifoldByStepOff(periodicOrbit, "Unstable", "Negative", 2*d, 10)
    @test isapprox(real.(UNManifold.initialConditions[2]), [0.8233687758541708, 4.0652881196013047E-11, 0, -3.4190188468627025E-10, 0.12655673501955023, 0])
    SPManifold::MBD.CR3BPManifold = getManifoldByTime(periodicOrbit, "Stable", "Positive", d, 10)
    @test (isapprox(real.(SPManifold.initialConditions[2]), [0.8530753931180122, 0.03396238898823798, 0, 0.013688866088412027, -0.10459166633216246, 0]) || isapprox(real.(SPManifold.initialConditions[2]), [0.8531910929863838, 0.034021565782213235, 0, 0.013326781331767417, -0.10470393601812733, 0]))
    @test isApproxSigFigs([getStabilityIndex(periodicOrbit)], [5.569779482796515E6], 8)
    @test isApproxSigFigs([getTimeConstant(periodicOrbit)], [0.35321154312125946], 8)
    @test getJacobiConstant(UPManifold) == getJacobiConstant(periodicOrbit)
    @test getJacobiConstant(manifoldArc) == getJacobiConstant(periodicOrbit)
end

@testset "SPICE Functions" begin
    SPICE.furnsh("../src/spice/kernels/naif0012.tls", "../src/spice/kernels/de440.bsp", "../src/spice/kernels/mar097.bsp")
    (ephemerisStates, ephemerisTimes) = getEphemerides("Jan 1 2026", [0.0, 2.743], "Earth", "Sun", "ECLIPJ2000")
    @test isApproxSigFigs(ephemerisStates[2], [-2.6074281011547867e7, 1.447742856780263e8, -8892.893132835627, -29.788851904426515, -5.3967021850401045, 0.0004108878034403407], 11)
    @test isApproxSigFigs([ephemerisTimes[2]], [490.68471150562976], 11)
    SPICE.kclear()
end

@testset "Utility Functions" begin
    @test isApproxSigFigs(Cartesian2Cylindrical([0.8234, 0.125, 0.4], [0.01, 0, 0]), [0.8229486982795, 0.1524830343520, 0.4], 13)
    @test_throws BoundsError checkIndices([1, 7], 6)
    @test_throws ArgumentError checkIndices([1, 1], 6)
    checkIndices([1, 5], 6)
    @test_throws ArgumentError maskData([true], [1.0 2.0 3.0; 4.0 5.0 6.0])
    mask1::Vector{Bool} = [true, false, true]
    data::Matrix{Float64} = [1.0 2.0 3.0; 4.0 5.0 6.0]
    @test maskData(mask1, data) == [1.0 3.0; 4.0 6.0]
    mask2::Vector{Bool} = [false, false, false]
    @test size(maskData(mask2, data)) == (1, 0)
    mask3::Vector{Bool} = [true, true, true]
    @test maskData(mask3, data) == data
end

# @testset "TBPDynamicsModel" begin
#     systemData = MBD.TBPSystemData("Earth")
#     dynamicsModel = MBD.TBPDynamicsModel(systemData)
#     @test_throws ArgumentError appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1.0], MBD.SIMPLE)
#     @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.SIMPLE) == [0.8234, 0, 0, 0, 0.1263, 0]
#     @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL) == [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
#     @test appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.ARCLENGTH) == [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0]
#     @test isApproxSigFigs(evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.8234, 0, 0, 0, 0.1263, 0]), [0, 0.1263, 0, -1.474953316253, 0, 0], 13)
#     @test isApproxSigFigs(evaluateEquations(dynamicsModel, MBD.FULL, 0.0, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL)), [0, 0.1263, 0, -1.474953316253, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0, 3.582592461143, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0], 13)
#     @test isApproxSigFigs(evaluateEquations(dynamicsModel, MBD.ARCLENGTH, 0.0, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.ARCLENGTH)), [0, 0.1263, 0, -1.474953316253, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0, 3.582592461143, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0, 0.1263], 13)
#     stateTrajectory::MBD.TBPTrajectory = getCartesianState(dynamicsModel, 100000.0, 0.1, 0.05, 1.0*pi, 1.0*pi, 0.0)
#     @test isApproxSigFigs(stateTrajectory.initialCondition, [90000.0, -2.202986797819E-11, 5.508614670424E-13, 5.402731156714E-16, 2.204453172147, -0.1103146027681], 13)
#     @test_throws ArgumentError getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
#     @test getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [0, 0, 0, 0, 0, 0]
#     @test isApproxSigFigs([getExcursion(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])], [1.231789044173E8], 13)
#     @test isApproxSigFigs(getLambertArc(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], 2.734, "Short")[1], [962.62693, -0.81942927, 0, 7.2123647, 144.41729, 0], 8)
#     @test isApproxSigFigs(getLambertArc(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], 2.734, "Short")[2], [-958.88199, 2.4108220, 0, -21.219315, -137.55345, 0], 8)
#     @test isApproxSigFigs(getLambertArc(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], 2.734, "Long")[1], [-313.30686, 96.924284, 0, -853.09776, 335.00117, 0], 8)
#     @test isApproxSigFigs(getLambertArc(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], 2.734, "Long")[2], [-297.06180, 96.896974, 0, -852.85739, 337.38504, 0], 8)
#     elementTrajectory::MBD.TBPTrajectory = getOsculatingOrbitalElements(dynamicsModel, stateTrajectory.initialCondition)
#     @test isApproxSigFigs([elementTrajectory.a], [100000.0], 13)
#     @test_throws ArgumentError getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
#     @test size(getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])) == (6, 0)
#     @test isApproxSigFigs([getPeriod(dynamicsModel, stateTrajectory)], [314710.3195678], 13)
#     @test getPrimaryPosition(dynamicsModel) == [0.0, 0.0, 0.0]
#     Moon = MBD.BodyData("Moon")
#     (resonantOrbit::MBD.TBPTrajectory, period::Float64) = getResonantOrbit(dynamicsModel, Moon, 3, 4, 0.7113)
#     @test isApproxSigFigs(resonantOrbit.initialCondition, [134559.8941719, 0, 0, 0, 2.251511352682, 0], 13)
#     @test getStateSize(dynamicsModel, MBD.SIMPLE) == 6
#     @test getStateSize(dynamicsModel, MBD.FULL) == 42
#     @test getStateSize(dynamicsModel, MBD.ARCLENGTH) == 43
#     @test_throws ArgumentError getStateTransitionMatrix(dynamicsModel, [0.8324, 0, 0, 0, 0.1263, 0])
#     @test getStateTransitionMatrix(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0]) == [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
#     @test isEpochIndependent(dynamicsModel)
#     @test isApproxSigFigs(primaryInertial2Rotating(dynamicsModel, Moon, [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[1], [0.8234, 0, 0, 0, 0.1262978217125, 0], 13)
#     @test isApproxSigFigs(primaryInertial2Rotating(dynamicsModel, Moon, [[0.8234, 0, 0, 0, 0.1263, 0], [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]], [0, 2.743])[2], [0.8322038162986, -0.002802017407416, 0, 0.02461018219822, 0.1165979143175, 0], 13)
#     KeplerTrajectory::MBD.TBPTrajectory = getCartesianState(dynamicsModel, 100000.0, 0.1, 0.05, 1.0*pi, 1.0*pi, 1.5)
#     @test isApproxSigFigs([solveKeplersEquation(dynamicsModel, KeplerTrajectory)], [65208.33705602], 13)
# end

# @testset "TBPEquationsOfMotion" begin
#     systemData = MBD.TBPSystemData("Earth")
#     dynamicsModel = MBD.TBPDynamicsModel(systemData)
#     EOMs_simple = MBD.TBPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
#     qdot_simple = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.SIMPLE))
#     computeDerivatives!(qdot_simple, [0.8234, 0, 0, 0, 0.1263, 0], (EOMs_simple,), 0.0)
#     @test isApproxSigFigs(qdot_simple, [0, 0.1263, 0, -1.474953316253, 0, 0], 13)
#     EOMs_full = MBD.TBPEquationsOfMotion(MBD.FULL, dynamicsModel)
#     qdot_full = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.FULL))
#     computeDerivatives!(qdot_full, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL), (EOMs_full,), 0.0)
#     @test isApproxSigFigs(qdot_full, [0, 0.1263, 0, -1.474953316253, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 3.582592461143, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0], 13)
#     EOMs_arclength = MBD.TBPEquationsOfMotion(MBD.ARCLENGTH, dynamicsModel)
#     qdot_arclength = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.ARCLENGTH))
#     computeDerivatives!(qdot_arclength, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.ARCLENGTH), (EOMs_arclength,), 0.0)
#     @test isApproxSigFigs(qdot_arclength, [0, 0.1263, 0, -1.474953316253, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 3.582592461143, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0, 0, 0, 0, -1.791296230572, 0, 0, 0, 0.1263], 13)
# end

# @testset "TBPTrajectory" begin
#     systemData = MBD.TBPSystemData("Earth")
#     dynamicsModel = MBD.TBPDynamicsModel(systemData)
#     stateTrajectory::MBD.TBPTrajectory = getCartesianState(dynamicsModel, 100000.0, 0.1, 0.05, 1.0*pi, 1.0*pi, 0.0)
#     KeplerTrajectory::MBD.TBPTrajectory = getCartesianState(dynamicsModel, 100000.0, 0.1, 0.05, 1.0*pi, 1.0*pi, 1.5)
#     @test getCartesianState(stateTrajectory, 1.5) == KeplerTrajectory
# end

@testset "Periodic Orbit Family Example" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    Moon = MBD.BodyData("Moon")
    L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
    (stateGuess::Vector{Float64}, timeGuess::Vector{Float64}) = getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])
    targetJC::Float64 = getJacobiConstant(dynamicsModel, stateGuess)
    targeter = LyapunovJCTargeter(dynamicsModel)
    solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuess, [timeGuess[1], timeGuess[2]], targetJC)
    @test LinearAlgebra.norm(getConstraintVector!(solution1)) <= 1E-11
    period1::Float64 = getPeriod(targeter, solution1)
    monodromy1::Matrix{Float64} = getMonodromy(targeter, solution1)
    orbit1 = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution1.nodes[1].state.data[1:6], period1, monodromy1)
    solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuess, [timeGuess[1], timeGuess[2]], targetJC-0.0001)
    @test LinearAlgebra.norm(getConstraintVector!(solution2)) <= 1E-11
    period2::Float64 = getPeriod(targeter, solution2)
    monodromy2::Matrix{Float64} = getMonodromy(targeter, solution2)
    orbit2 = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution1.nodes[1].state.data[1:6], period2, monodromy2)
    JCContinuation = MBD.JacobiConstantContinuationEngine(solution1, solution2, -1E-3, -1E-2)
    JCContinuation.printProgress = false
    ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [NaN NaN; -2.5 0])
    addJumpCheck!(JCContinuation, ydot0JumpCheck)
    stepsEndCheck = MBD.NumberStepsContinuationEndCheck(100)
    MoonEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [L1[1] getPrimaryPosition(dynamicsModel, 2)[1]-Moon.bodyRadius/getCharLength(systemData); NaN NaN])
    addEndCheck!(JCContinuation, stepsEndCheck)
    addEndCheck!(JCContinuation, MoonEndCheck)
    oldstd = stdout
    redirect_stdout(devnull)
    solutions::MBD.CR3BPContinuationFamily = doContinuation!(JCContinuation, solution1, solution2)
    @test getNumMembers(solutions) == 100
    orbit100::MBD.CR3BPPeriodicOrbit = getPeriodicOrbit(targeter, solutions, 100)
    @test isApproxSigFigs([getJacobiConstant(orbit100)], [2.9907766923863504], 8)
    orbitFamily = MBD.CR3BPOrbitFamily(dynamicsModel)
    orbitFamily2 = MBD.CR3BPOrbitFamily(dynamicsModel)
    @test orbitFamily == orbitFamily2
    for s::Int64 in 1:getNumMembers(solutions)
        orbit::MBD.CR3BPPeriodicOrbit = getPeriodicOrbit(targeter, solutions, s)
        push!(orbitFamily.initialConditions, orbit.initialCondition)
        push!(orbitFamily.periods, orbit.period)
        push!(orbitFamily.monodromies, orbit.monodromy)
    end
    @test getNumMembers(orbitFamily) == getNumMembers(solutions)
    @test getJacobiConstant(dynamicsModel, orbitFamily.initialConditions[100]) == getJacobiConstant(orbit100)
    redirect_stdout(devnull)
    eigenSort!(orbitFamily)
    redirect_stdout(oldstd)
    @test isApproxSigFigs(real.(orbitFamily.eigenvalues[100]), [0.004034714876274129, 247.8489882397049, 0.5901934719701916, 0.5901934719701916, 0.999994007551255, 1.0000059924982916], 4)
    (standardIndices::Vector{Vector{Float64}}, alternateIndices::Vector{Vector{Complex{Float64}}}) = getAlternateIndices(orbitFamily)
    @test isApproxSigFigs(standardIndices[100], [123.9265115235783, 1.0, 1.0000000000282148], 5)
    @test isApproxSigFigs(real.(alternateIndices[100]), [123.92651152292481, 0.5901934719700709, 1.0000000000255436], 5)
    x0Targeter = Lyapunovx0Targeter(dynamicsModel)
    solution3::MBD.CR3BPMultipleShooterProblem = correct(x0Targeter, stateGuess, [timeGuess[1], timeGuess[2]], stateGuess[1])
    solution4::MBD.CR3BPMultipleShooterProblem = correct(x0Targeter, stateGuess, [timeGuess[1], timeGuess[2]], stateGuess[1]+1E-4)
    NPContinuation = MBD.CR3BPNaturalParameterContinuationEngine(solution3, solution4, "Initial State", 1, 1E-4, 1E-2)
    NPContinuation.printProgress = false
    addJumpCheck!(NPContinuation, ydot0JumpCheck)
    addEndCheck!(NPContinuation, stepsEndCheck)
    addEndCheck!(NPContinuation, MoonEndCheck)
    redirect_stdout(devnull)
    solutions2::MBD.CR3BPContinuationFamily = doContinuation!(NPContinuation, solution3, solution4)
    redirect_stdout(oldstd)
    @test getNumMembers(solutions2) == 100
    newOrbit100::MBD.CR3BPPeriodicOrbit = getPeriodicOrbit(x0Targeter, solutions2, 100)
    @test isApproxSigFigs([newOrbit100.initialCondition[1]], [0.8616151323643001], 8)
#     detectBifurcations!(targeter, family)
#     bifurcation::MBD.Bifurcation = family.bifurcations[1]
#     @test isApproxSigFigs(bifurcation.orbit.initialCondition, [0.85479945002, 0, 0, 0, -0.13373284140, 0], 11)
end

@testset "Earth-Moon L1 Halo Multiple Shooter Example" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    nodeStates::Vector{Vector{Float64}} = [[1.0277, 0, 0.1857, 0, -0.1152, 0], [1.008, -0.0476, 0.1133, -0.0714, -0.0307, -0.2984], [1.008, 0.0476, 0.1133, 0.0714, -0.0307, 0.2985], [1.0277, 0, 0.1857, 0, -0.1152, 0]]
    nodeTimes::Vector{Float64} = [0, 0.52854, 1.0571, 1.5856]
    nodes::Vector{MBD.CR3BPNode} = Vector{MBD.CR3BPNode}(undef, 4)
    [nodes[n] = MBD.CR3BPNode(nodeTimes[n], nodeStates[n], dynamicsModel) for n = 1:4]
    problem = MBD.CR3BPMultipleShooterProblem()
    segments::Vector{MBD.CR3BPSegment} = Vector{MBD.CR3BPSegment}(undef, 3)
    for s::Int16 in Int16(1):Int16(3)
        segments[s] = MBD.CR3BPSegment(nodeTimes[s+1]-nodeTimes[s], nodes[s], nodes[s+1])
        addSegment!(problem, segments[s])
    end
    map(s -> addConstraint!(problem, MBD.CR3BPContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.StateMatchConstraint(nodes[1].state, nodes[end].state, [1, 2, 3, 4, 5, 6]))
    addConstraint!(problem, MBD.JacobiConstraint(nodes[1], 3.04))
    @test checkJacobian(problem)
    multipleShooter = MBD.CR3BPMultipleShooter()
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(multipleShooter, problem)
    @test isApproxSigFigs(solution.nodes[1].state.data, solution.nodes[end].state.data, 11)
end
