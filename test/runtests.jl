using MBD
using Test

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
include("../src/propagation/Propagator.jl")
include("../src/spice/BodyName.jl")
include("../src/utilities/UtilityFunctions.jl")

@testset "Files" begin
    @test isfile("../src/body_data.xml")
end

@testset "Constructors" begin
    @test MBD.CR3BPSystemData <: MBD.AbstractSystemData
    @test_throws ArgumentError MBD.CR3BPSystemData("Earth", "Mars")
    @test MBD.CR3BPDynamicsModel <: MBD.AbstractDynamicsModel
    @test MBD.CR3BPEquationsOfMotion <: MBD.AbstractEquationsOfMotion
    @test_throws ArgumentError MBD.Variable([1.0, 1.0], [true, false, true])
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    @test MBD.Node <: MBD.IHasVariables
    @test_throws ArgumentError MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263], dynamicsModel)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    @test MBD.Segment <: MBD.IHasVariables
    @test_throws ArgumentError MBD.Segment(2.743, originNode, terminalNode)
    @test MBD.MultipleShooterProblem <: MBD.AbstractNonlinearProblem
    @test MBD.ContinuityConstraint <: MBD.AbstractConstraint
    @test MBD.StateConstraint <: MBD.AbstractConstraint
    @test MBD.StateMatchConstraint <: MBD.AbstractConstraint
    variable1 = MBD.Variable([1.0, 2.0], [true, false])
    variable2 = MBD.Variable([1.0, 2.0, 3.0], [true, false, true])
    @test_throws ArgumentError MBD.StateMatchConstraint(variable1, variable2, [1, 2])
    variable3 = MBD.Variable([1.0, 2.0], [true, true])
    @test_throws ArgumentError MBD.StateMatchConstraint(variable1, variable3, [1, 2])
    @test MBD.JacobiConstraint <: MBD.AbstractConstraint
    @test MBD.ConstraintVectorL2NormConvergenceCheck <: MBD.AbstractConvergenceCheck
    @test MBD.MinimumNormUpdateGenerator <: MBD.AbstractUpdateGenerator
    @test MBD.LeastSquaresUpdateGenerator <: MBD.AbstractUpdateGenerator
    @test MBD.MultipleShooter <: MBD.AbstractNonlinearProblemSolver
    @test_throws ArgumentError MBD.AdaptiveStepSizeByElementGenerator("Initial State", -1, 1E-4, 1E-2)
    @test MBD.NumberStepsContinuationEndCheck <: MBD.AbstractContinuationEndCheck
    @test MBD.BoundingBoxContinuationEndCheck <: MBD.AbstractContinuationEndCheck
    @test MBD.BoundingBoxJumpCheck <: MBD.AbstractContinuationJumpCheck
    @test MBD.NaturalParameterContinuationEngine <: MBD.AbstractContinuationEngine
    @test MBD.JacobiConstantContinuationEngine <: MBD.AbstractContinuationEngine
    @test MBD.CR3BPPeriodicOrbit <: MBD.AbstractTrajectoryStructure
end

@testset "Copy" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    @test MBD.shallowClone(originNode) == originNode
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    @test MBD.shallowClone(segment) == segment
    continuityConstraint = MBD.ContinuityConstraint(segment)
    @test MBD.shallowClone(continuityConstraint) == continuityConstraint
    stateConstraint = MBD.StateConstraint(originNode, [1], [originNode.state.data[1]])
    @test MBD.shallowClone(stateConstraint) == stateConstraint
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test MBD.shallowClone(stateMatchConstraint) == stateMatchConstraint
    jacobiConstraint = MBD.JacobiConstraint(originNode, 3.1743560232059265)
    @test MBD.shallowClone(jacobiConstraint) == jacobiConstraint
end

@testset "Deep Copy" begin
    variable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, false, false])
    @test MBD.deepClone(variable) == variable
end

@testset "BodyName" begin
    bodyName = MBD.BodyName("Earth")
    @test getIDCode(bodyName) == 399
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
    @test getLinearVariation(dynamicsModel, L1, [0.005, 0, 0]) == ([0.8419151323643023, 0, 0, 0, -0.04186136597442648, 0], [0, 2.6915795607981865])
    @test_throws ArgumentError getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test size(getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])) == (6, 0)
    @test_throws ArgumentError getPrimaryPosition(dynamicsModel, 3)
    @test getPrimaryPosition(dynamicsModel, 1) == [-0.012150584269940356, 0, 0]
    @test getPseudopotentialJacobian(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0]) == [9.851145859594304, -3.4255729297971507, -4.42557292979715, 0, 0, 0]
    @test_throws ArgumentError getStateTransitionMatrix(dynamicsModel, [0.8324, 0, 0, 0, 0.1263, 0])
    @test getStateTransitionMatrix(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
    @test_throws ArgumentError primaryInertial2Rotating(dynamicsModel, 2, [[1.0, 0, 0, 0, 1.0, 0], [1.1, 0, 0, 0, 1.1, 0]], [0.0])
    @test_throws ArgumentError primaryInertial2Rotating(dynamicsModel, 3, [[1.0, 0, 0, 0, 1.0, 0]], [0.0])
end

@testset "CR3BPEquationsOfMotion" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    EOMs_simple = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    qdot_simple = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.SIMPLE))
    computeDerivatives!(qdot_simple, [0.8234, 0, 0, 0, 0.1263, 0], EOMs_simple, 0.0)
    @test qdot_simple == [0, 0.1263, 0, 0.11033238649399063, 0, 0]
    EOMs_full = MBD.CR3BPEquationsOfMotion(MBD.FULL, dynamicsModel)
    qdot_full = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.FULL))
    computeDerivatives!(qdot_full, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL), EOMs_full, 0.0)
    @test qdot_full == [0, 0.1263, 0, 0.11033238649399063, 0, 0, 0, 0, 0, 9.851145859594295, 0, 0, 0, 0, 0, 0, -3.425572929797148, 0, 0, 0, 0, 0, 0, -4.4255729297971484, 1, 0, 0, 0, -2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0]
end

@testset "Propagator" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.Arc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test arc.states[end] == [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]
    @test arc.times[end-1] == 2.6982233859741345
end

@testset "Arc" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.Arc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test getStateByIndex(arc, -1) == [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]
    @test getTimeByIndex(arc, -1) == 2.743
    deleteStateAndTime!(arc, 45)
    @test getStateCount(arc) == 45
    @test arc.states[45] == [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]
    @test arc.times[45] == 2.743
    @test_throws BoundsError deleteStateAndTime!(arc, 46)
    @test_throws BoundsError getStateByIndex(arc, 46)
    @test_throws BoundsError getTimeByIndex(arc, 46)
    setParameters!(arc, systemData.params)
    @test arc.params == systemData.params
end

@testset "Variable" begin
    variable = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, false, false, false, true, false])
    @test getData(variable) == [0.8234, 0, 0, 0, 0.1263, 0]
    @test getFreeVariableMask(variable) == [true, false, false, false, true, false]
    @test getFreeVariableData(variable) == [0.8234, 0.1263]
    @test getNumberFreeVariables(variable) == 2
    @test_throws ArgumentError setFreeVariableData!(variable, [4.0, 7.0, 8.5])
    setFreeVariableData!(variable, [4.0, 7.0])
    @test variable.data == [4.0, 0, 0, 0, 7.0, 0]
    @test_throws ArgumentError setFreeVariableMask!(variable, [true, false])
    setFreeVariableMask!(variable, [true, false, false, false, false, false])
    @test variable.freeVariableMask == [true, false, false, false, false, false]
end

@testset "Node" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    node = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    state = MBD.Variable([0.8234, 0, 0, 0, 0.1263, 0], [true, true, true, true, true, true])
    state.name = "Node State"
    epoch = MBD.Variable([0.0], [false])
    epoch.name = "Node Epoch"
    @test getVariables(node) == [state, epoch]
end

@testset "Segment" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    @test getPropagatorParametersData(segment) == []
    lazyPropagate!(segment, MBD.SIMPLE)
    @test getFinalState!(segment) == [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]
    @test getFinalStateRate!(segment) == [0.024609343495764855, 0.1166002944773056, 0, 0.18113477239159148, -0.03842090365548907, 0]
    @test getPartials_FinalStateWRTEpoch!(segment) == zeros(Float64, (6, 1))
    @test getPartials_FinalStateWRTInitialState!(segment) == [1363.7623463462148 -350.51893299536675 0 401.47006684424275 130.91187990975433 0; -421.081159205987 109.07625516626688 0 -123.8052684899764 -40.551373197002775 0; 0 0 0.9853225289729513 0 0 -0.07418084339056366; 3877.9818867917584 -995.8827212743846 0 1141.6513906028536 372.4229726062542 0; -1465.2146870834852 376.8074206167876 0 -431.6882797615755 -139.73106841637286 0; 0 0 -0.09704680345808286 0 0 1.022202359240349]
    @test size(getPartials_FinalStateWRTParams!(segment)) == (6, 0)
    lazyPropagate!(segment, MBD.FULL)
    @test getFinalState!(segment) == [0.8322038366146268, -0.0027959784808503307, 0, 0.02460934350980587, 0.11660029447199295, 0, 1363.7623463462148, -421.081159205987, 0, 3877.9818867917584, -1465.2146870834852, 0, -350.51893299536675, 109.07625516626688, 0, -995.8827212743846, 376.8074206167876, 0, 0, 0, 0.9853225289729513, 0, 0, -0.09704680345808286, 401.47006684424275, -123.8052684899764, 0, 1141.6513906028536, -431.6882797615755, 0, 130.91187990975433, -40.551373197002775, 0, 372.4229726062542, -139.73106841637286, 0, 0, 0, -0.07418084339056366, 0, 0, 1.022202359240349]
    @test getVariables(segment) == [segment.TOF, segment.propParams]
    resetPropagatedArc!(segment)
    @test getStateCount(segment.propArc) == 1
end

@testset "ContinuityConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    segment = MBD.Segment(2.743, originNode, terminalNode)
    continuityConstraint = MBD.ContinuityConstraint(segment)
    @test getNumberConstraintRows(continuityConstraint) == 6
    multipleShooterProblem = MBD.MultipleShooterProblem()
    @test evaluateConstraint(continuityConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector) == [0, 0, 0, 0, 0, 0]
    @test length(keys(getPartials_ConstraintWRTVariables(continuityConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector))) == 3
end

@testset "StateConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    stateConstraint = MBD.StateConstraint(originNode, [1], [originNode.state.data[1]])
    @test getNumberConstraintRows(stateConstraint) == 1
    multipleShooterProblem = MBD.MultipleShooterProblem()
    @test evaluateConstraint(stateConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector) == [0]
    @test length(keys(getPartials_ConstraintWRTVariables(stateConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector))) == 1
end

@testset "StateMatchConstraint" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    originNode = MBD.Node(0.0, [0.8234, 0, 0, 0, 0.1263, 0], dynamicsModel)
    terminalNode = MBD.Node(2.743, [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0], dynamicsModel)
    stateMatchConstraint = MBD.StateMatchConstraint(originNode.state, terminalNode.state, [1, 2, 3, 4, 5, 6])
    @test getNumberConstraintRows(stateMatchConstraint) == 6
    multipleShooterProblem = MBD.MultipleShooterProblem()
    @test evaluateConstraint(stateMatchConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector) == [-0.008803836609691418, 0.00279597847933625, 0, -0.024609343495764855, 0.00969970552269439, 0]
    @test length(keys(getPartials_ConstraintWRTVariables(stateMatchConstraint, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector))) == 2
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
    @test qf == [0.8322038366096914, -0.00279597847933625, 0, 0.024609343495764855, 0.1166002944773056, 0]
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
