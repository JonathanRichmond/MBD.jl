using MBD
using Test

include("../src/spice/BodyName.jl")

@testset "Constructors" begin
    @test MBD.SystemData <: MBD.AbstractSystemData
    @test_throws ArgumentError MBD.SystemData("Earth", "Mars")
end
