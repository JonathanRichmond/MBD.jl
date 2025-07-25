"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 6/25/25
"""

using Logging, MBD, Test
using DifferentialEquations, LightXML, SPICE

global_logger(ConsoleLogger(stderr, Logging.Info)) # Debug/Info/Warn/Error

@testset "Files" begin
    @test isfile("../src/body_data.xml")
    # @test isfile("../src/spice/kernels/de440.bsp")
    # @test isfile("../src/spice/kernels/naif0012.tls")
end

@testset "BodyData Tests" begin
    @testset "Constructor with valid body name" begin
        bodyData = MBD.BodyData("Earth")
        println(bodyData)
        display(bodyData)

        @test typeof(bodyData) == MBD.BodyData
        @test bodyData.name == "Earth"
        @test bodyData.SPICEID isa Int16
        @test bodyData.gravParam > 0.0
        @test bodyData.mass > 0.0
        @test bodyData.bodyRadius > 0.0
        @test bodyData.orbitRadius > bodyData.bodyRadius
        @test bodyData.inc >= 0.0
        @test bodyData.RAAN >= 0.0
    end

    @testset "Invalid SPICE name" begin
        err = try
            MBD.BodyData("FakePlanet")
            nothing
        catch e
            e
        end

        @test err isa MethodError
        @test occursin("no method matching", sprint(showerror, err))
    end

    @testset "Valid SPICE ID not in XML" begin
        fakeName::String = "Daphnis"
        if try
            SPICE.bods2c(fakeName)
            true
        catch
            false
        end
        err = try
            MBD.BodyData(fakeName)
            nothing
        catch e
            e
        end

        @test err isa ArgumentError
        @test occursin("not found in", sprint(showerror, err))
        else
            Logging.@info "Skipping test: SPICE ID for '$fakeName' not defined in current kernel set"
        end
    end

    @testset "Corrupted XML" begin
        xml_path::String = joinpath(@__DIR__, "..", "src", "body_data.xml")
        corrupt_path::String = xml_path*".bak"
        try
            isfile(xml_path) || Logging.@warn "Test skipped: body_data.xml not found"
            mv(xml_path, corrupt_path)
            err = try
                MBD.BodyData("Earth")
                nothing
            catch e
                e
            end

            @test err isa LightXML.XMLParseError
            @test occursin("Failure in parsing", sprint(showerror, err))
        finally
            isfile(corrupt_path) && mv(corrupt_path, xml_path)
        end
    end

    @testset "Equality" begin
        bodyData1 = MBD.BodyData("Earth")
        bodyData2 = MBD.BodyData("Earth")
        bodyData3 = MBD.BodyData("Moon")

        @test bodyData1 == bodyData2
        @test bodyData1 != bodyData3
    end
end

@testset "SystemData Tests" begin
    @testset "Constructor with valid 3-body system" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        println(systemData)
        display(systemData)

        @test typeof(systemData) == MBD.SystemData
        @test systemData.modelType == MBD.CR3BP
        @test length(systemData.primaryNames) == 2
        @test systemData.primaryNames == ["Earth", "Moon"]
        @test all(x -> x isa MBD.BodyData, systemData.primaryData)
        @test all(x -> x isa Int16, systemData.primarySPICEIDs)
    end

    @testset "Constructor with valid 4-body system" begin
        systemData = MBD.SystemData(MBD.BCR4BP, "Sun", "Earth", "Moon")

        @test typeof(systemData) == MBD.SystemData
        @test systemData.modelType == MBD.BCR4BP
        @test length(systemData.primaryNames) == 3
        @test systemData.primaryNames == ["Sun", "Earth", "Moon"]
        @test all(x -> x isa MBD.BodyData, systemData.primaryData)
        @test all(x -> x isa Int16, systemData.primarySPICEIDs)
    end

    @testset "Invalid model type" begin
        err = try
            MBD.SystemData(MBD.QBCR4BP, "Sun", "Earth", "Moon")
            nothing
        catch e
            e
        end

        @test err isa UndefVarError
        @test occursin("not defined in", sprint(showerror, err))
    end

    @testset "Invalid body name" begin
        err = try
            MBD.SystemData(MBD.CR3BP, "Earth", "FakePlanet")
            nothing
        catch e
            e
        end

        @test err isa Exception
        @test occursin("no method matching", sprint(showerror, err)) || occursin("not found in", sprint(showerror, err)) || occursin("Failure in parsing", sprint(showerror, err))
    end

    @testset "Improper parent relationship" begin
        @test_logs match_mode = :any min_level = Logging.Info (
            (:warn, r".* does not have .* as parent")
        ) begin
            MBD.SystemData(MBD.CR3BP, "Jupiter", "Earth")
        end
    end

    @testset "Empty system" begin
        systemData = MBD.SystemData(MBD.TBP)

        @test typeof(systemData) == MBD.SystemData
        @test systemData.modelType == MBD.TBP
        @test isempty(systemData.primaryNames)
        @test isempty(systemData.primaryData)
        @test isempty(systemData.primarySPICEIDs)
    end

    @testset "Equality" begin
        systemData1 = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        systemData2 = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        systemData3 = MBD.SystemData(MBD.CR3BP, "Sun", "Earth")

        @test systemData1 == systemData2
        @test systemData1 != systemData3
    end
end

@testset "DynamicsModel Tests" begin
    @testset "Constructor with valid system data" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        dynamicsModel = MBD.DynamicsModel(systemData)
        println(dynamicsModel)
        display(dynamicsModel)

        @test typeof(dynamicsModel) == MBD.DynamicsModel
        @test dynamicsModel.systemData === systemData
    end

    @testset "Equality" begin
        systemData1 = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        systemData2 = MBD.SystemData(MBD.CR3BP, "Sun", "Earth")
        dynamicsModel1 = MBD.DynamicsModel(systemData1)
        dynamicsModel2 = MBD.DynamicsModel(systemData1)
        dynamicsModel3 = MBD.DynamicsModel(systemData2)

        @test dynamicsModel1 == dynamicsModel2
        @test dynamicsModel1 != dynamicsModel3
    end
end

@testset "EquationsOfMotion Tests" begin
    @testset "Constructor with valid simple dynamics model" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        dynamicsModel = MBD.DynamicsModel(systemData)
        equations = MBD.EquationsOfMotion(dynamicsModel, MBD.SIMPLE)
        println(equations)
        display(equations)

        @test typeof(equations) == MBD.EquationsOfMotion
        @test equations.equationType == MBD.SIMPLE
        @test equations.dynamicsModel === dynamicsModel
    end

    @testset "Constructor with valid STM dynamics model" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        dynamicsModel = MBD.DynamicsModel(systemData)
        equations = MBD.EquationsOfMotion(dynamicsModel, MBD.STM)

        @test typeof(equations) == MBD.EquationsOfMotion
        @test equations.equationType == MBD.STM
    end

    @testset "Invalid equation type" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        dynamicsModel = MBD.DynamicsModel(systemData)
        err = try
            MBD.EquationsOfMotion(dynamicsModel, MBD.COMPLEX)
            nothing
        catch e
            e
        end

        @test err isa UndefVarError
        @test occursin("not defined in", sprint(showerror, err))
    end

    @testset "Equality" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        dynamicsModel = MBD.DynamicsModel(systemData)
        equations1 = MBD.EquationsOfMotion(dynamicsModel, MBD.SIMPLE)
        equations2 = MBD.EquationsOfMotion(dynamicsModel, MBD.SIMPLE)
        equations3 = MBD.EquationsOfMotion(dynamicsModel, MBD.STM)

        @test equations1 == equations2
        @test equations1 != equations3
    end
end

@testset "Integrator Tests" begin
    @testset "Constructor with valid VERN9 integrator" begin
        integrator = MBD.Integrator(MBD.VERN9)
        println(integrator)
        display(integrator)

        @test typeof(integrator) == MBD.Integrator
        @test integrator.integratorType == MBD.VERN9
    end

    @testset "Constructor with valid DP8 integrator" begin
        integrator = MBD.Integrator(MBD.DP8)

        @test typeof(integrator) == MBD.Integrator
        @test integrator.integratorType == MBD.DP8
    end

    @testset "Invalid integrator type" begin
        err = try
            MBD.Integrator(MBD.VERN11)
            nothing
        catch e
            e
        end

        @test err isa UndefVarError
        @test occursin("not defined in", sprint(showerror, err))
    end

    @testset "Equality" begin
        integrator1 = MBD.Integrator(MBD.VERN9)
        integrator2 = MBD.Integrator(MBD.VERN9)
        integrator3 = MBD.Integrator(MBD.DP8)

        @test integrator1 == integrator2
        @test integrator1 != integrator3
    end
end

@testset "Propagator Tests" begin
    @testset "Constructor with defaults" begin
        propagator = MBD.Propagator()
        println(propagator)
        display(propagator)

        @test typeof(propagator) == MBD.Propagator
        @test propagator.absTol == 1E-12
        @test propagator.relTol == 1E-12
        @test propagator.maxStep == 100
        @test propagator.maxEvaluations == typemax(Int64)
        @test propagator.events == []
        @test propagator.equationType == MBD.SIMPLE
        @test propagator.integrator.integratorType == MBD.DP8
    end

    @testset "Constructor with valid VERN9 integrator and STM equation type" begin
        integrator = MBD.Integrator(MBD.VERN9)
        propagator = MBD.Propagator(integrator = integrator, equationType = MBD.STM)

        @test typeof(propagator) == MBD.Propagator
        @test propagator.integrator.integratorType == MBD.VERN9
        @test propagator.equationType == MBD.STM
    end

    @testset "Invalid equation type" begin
        err = try
            MBD.Propagator(equationType = MBD.COMPLEX)
            nothing
        catch e
            e
        end

        @test err isa UndefVarError
        @test occursin("not defined in", sprint(showerror, err))
    end

    @testset "Equality" begin
        propagator1 = MBD.Propagator()
        propagator2 = MBD.Propagator()

        @test propagator1 == propagator2
    end
end

@testset "Arc Tests" begin
    @testset "Constructor with valid dynamics model" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        dynamicsModel = MBD.DynamicsModel(systemData)
        arc = MBD.Arc(dynamicsModel)
        println(arc)
        display(arc)

        @test typeof(arc) == MBD.Arc
        @test arc.dynamicsModel === dynamicsModel
        @test arc.states == []
        @test arc.times == []
    end

    @testset "Equality" begin
        systemData1 = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")
        systemData2 = MBD.SystemData(MBD.CR3BP, "Sun", "Earth")
        dynamicsModel1 = MBD.DynamicsModel(systemData1)
        dynamicsModel2 = MBD.DynamicsModel(systemData2)
        arc1 = MBD.Arc(dynamicsModel1)
        arc2 = MBD.Arc(dynamicsModel1)
        arc3 = MBD.Arc(dynamicsModel2)

        @test arc1 == arc2
        @test arc1 != arc3
    end
end
