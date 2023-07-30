using MBD
using Test

include("../src/spice/BodyName.jl")

@testset "Files" begin
    @test isfile("../src/body_data.xml")
end

@testset "Constructors" begin
    @test MBD.CR3BPSystemData <: MBD.AbstractSystemData
    @test_throws ArgumentError MBD.CR3BPSystemData("Earth", "Mars")
    @test MBD.CR3BPDynamicsModel <: MBD.AbstractDynamicsModel
end

@testset "BodyName" begin
    @test getIDCode(MBD.BodyName("Earth")) == 399
end