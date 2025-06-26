"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 6/25/25
"""

using Logging, MBD, Test
using LightXML, SPICE

global_logger(ConsoleLogger(stderr, Logging.Info)) # Debug/Info/Warn/Error

@testset "Files" begin
    @test isfile("../src/body_data.xml")
    # @test isfile("../src/spice/kernels/de440.bsp")
    # @test isfile("../src/spice/kernels/naif0012.tls")
end

@testset "BodyData Tests" begin
    @testset "Constructor with valid body name" begin
        body = MBD.BodyData("Earth")

        @test typeof(body) == MBD.BodyData
        @test body.name == "Earth"
        @test body.SPICEID isa Int16
        @test body.gravParam > 0.0
        @test body.mass > 0.0
        @test body.bodyRadius > 0.0
        @test body.orbitRadius > body.bodyRadius
        @test body.inc >= 0.0
        @test body.RAAN >= 0.0
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
        body1 = MBD.BodyData("Earth")
        body2 = MBD.BodyData("Earth")

        @test body1 == body2
    end
end

@testset "SystemData Tests" begin
    @testset "Constructor with valid 3-body system" begin
        systemData = MBD.SystemData(MBD.CR3BP, "Earth", "Moon")

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

        @test systemData1 == systemData2
    end
end
