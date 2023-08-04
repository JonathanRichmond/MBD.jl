using MBD
using Test

include("../src/CR3BP/DynamicsModel.jl")
include("../src/CR3BP/EquationsOfMotion.jl")
include("../src/CR3BP/SystemData.jl")
include("../src/propagation/Propagator.jl")
include("../src/spice/BodyName.jl")

@testset "Files" begin
    @test isfile("../src/body_data.xml")
end

@testset "Constructors" begin
    @test MBD.CR3BPSystemData <: MBD.AbstractSystemData
    @test_throws ArgumentError MBD.CR3BPSystemData("Earth", "Mars")
    @test MBD.CR3BPDynamicsModel <: MBD.AbstractDynamicsModel
    @test MBD.CR3BPEquationsOfMotion <: MBD.AbstractEquationsOfMotion
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
    @test evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.8234, 0, 0, 0, 0.1263, 0], systemData.params) == [0, 0.1263, 0, 0.11033238649399063, -1.6468, 0]
    @test getEquationsOfMotion(dynamicsModel, MBD.SIMPLE) == MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    @test isEpochIndependent(dynamicsModel)
    @test_throws ArgumentError getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test getEpochDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [0, 0, 0, 0, 0, 0]
    @test_throws ArgumentError getEquilibriumPoint(dynamicsModel, getMassRatio(systemData), 6)
    @test getEquilibriumPoint(dynamicsModel, getMassRatio(systemData), 1) == [0.8369151323643023, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, getMassRatio(systemData), 2) == [1.1556821602923404, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, getMassRatio(systemData), 3) == [-1.005062645252109, 0, 0]
    @test getEquilibriumPoint(dynamicsModel, getMassRatio(systemData), 4) == [0.48784941573005963, 0.8660254037844386, 0]
    @test getEquilibriumPoint(dynamicsModel, getMassRatio(systemData), 5) == [0.48784941573005963, -0.8660254037844386, 0]
    @test getJacobiConstant(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0]) == 3.1743560232059265
    @test_throws ArgumentError getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0])
    @test size(getParameterDependencies(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])) == (6, 0)
    @test_throws ArgumentError getPrimaryPosition(dynamicsModel, getMassRatio(systemData), 3)
    @test getPrimaryPosition(dynamicsModel, getMassRatio(systemData), 1) == [-0.012150584269940356, 0, 0]
    @test getPseudopotentialJacobian(dynamicsModel, getMassRatio(systemData), [0.8234, 0, 0, 0, 0.1263, 0]) == [9.851145859594304, -3.4255729297971507, -4.42557292979715, 0, 0, 0]
    @test_throws ArgumentError getStateTransitionMatrix(dynamicsModel, [0.8324, 0, 0, 0, 0.1263, 0])
    @test getStateTransitionMatrix(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]) == [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
end

@testset "CR3BPEquationsOfMotion" begin
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    EOMs_simple = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, dynamicsModel)
    qdot_simple = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.SIMPLE))
    computeDerivatives!(qdot_simple, [0.8234, 0, 0, 0, 0.1263, 0], EOMs_simple, 0.0)
    @test qdot_simple == [0, 0.1263, 0, 0.11033238649399063, -1.6468, 0]
    EOMs_full = MBD.CR3BPEquationsOfMotion(MBD.FULL, dynamicsModel)
    qdot_full = Vector{Float64}(undef, getStateSize(dynamicsModel, MBD.FULL))
    computeDerivatives!(qdot_full, appendExtraInitialConditions(dynamicsModel, [0.8234, 0, 0, 0, 0.1263, 0], MBD.FULL), EOMs_full, 0.0)
    @test qdot_full == [0, 0.1263, 0, 0.11033238649399063, -1.6468, 0, 0, 0, 0, 9.851145859594295, 0, 0, 0, 0, 0, 0, -3.425572929797148, 0, 0, 0, 0, 0, 0, -4.4255729297971484, 1, 0, 0, 0, -2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0]
end

@testset "Propagator" begin
    propagator = MBD.Propagator()
    systemData = MBD.CR3BPSystemData("Earth", "Moon")
    dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
    arc::MBD.Arc = propagate(propagator, [0.8234, 0, 0, 0, 0.1263, 0], [0, 2.743], dynamicsModel)
    @test arc.states[60] == [-2.6985427953764827, 1.6414206860460856, 0, 0.5136679523797875, 4.9939829755428615, 0]
    @test arc.times[59] == 2.7348634602445325
end
