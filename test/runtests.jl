"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 12/12/25
"""

using MBD, Test

import LightXML, Logging, SPICE


Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info)) # Debug/Info/Warn/Error


@testset "Constructors" begin
    @testset "BodyData constructors" begin
        # Create a temporary XML file containing a single test body with a known ID
        xml = """<?xml version="1.0"?>
        <bodies>
        <body>
            <id>9999</id>
            <circ_r>7000.0</circ_r>
            <ecc>0.001</ecc>
            <inc>0.1</inc>
            <parentId>NaN</parentId>
            <radius>6371.0</radius>
            <gm>398600.4418</gm>
            <raan>0.5</raan>
        </body>
        </bodies>
        """

        tmpfile = tempname() * ".xml"
        open(tmpfile, "w") do io
            write(io, xml)
        end

        # Temporarily override MBD.getIDCode to return the test ID (avoids SPICE dependency)
        orig_getIDCode = MBD._getIDCode_func[]
        MBD._getIDCode_func[] = (_) -> 9999

        try
            bd = MBD.load_bodyData("TestBody", tmpfile)
            @test isa(bd, MBD.BodyData)
            @test isapprox(bd.a, 7000.0; atol=1e-12)
            @test isapprox(bd.e, 0.001; atol=1e-12)
            @test isapprox(bd.i, 0.1; atol=1e-12)
            @test bd.parentSPICEID == Int16(MBD.UNINITIALIZED_INDEX)
            @test isapprox(bd.r, 6371.0; atol=1e-12)
            @test isapprox(bd.μ, 398600.4418; atol=1e-9)
            @test bd.SPICEID == Int16(9999)

            # Ensure the pretty-print shows the name
            s = sprint(show, bd)
            @test occursin("BodyData:", s) && occursin("TestBody", s)
        finally
            # Restore original getIDCode and remove temp file
            MBD._getIDCode_func[] = orig_getIDCode
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @testset "SystemData constructors" begin
        # Map names to the IDs present in src/body_data.xml
        orig_resolver = MBD._getIDCode_func[]
        MBD._getIDCode_func[] = name -> begin
            n = lowercase(strip(name))
            if n == "earth"
                return 399
            elseif n == "moon"
                return 301
            else
                return 0
            end
        end

        try
            sys = MBD.init_systemData(["Earth", "Moon"])
            @test isa(sys, MBD.SystemData)
            @test length(sys.bodyData) == 2
            @test sys.names == ["Earth", "Moon"]
            @test sys.SPICEIDs[1] == Int16(399)
            @test sys.SPICEIDs[2] == Int16(301)

            # Ensure the pretty-print shows the name
            s = sprint(show, sys)
            @test occursin("SystemData:", s) && occursin("Earth", s)
        finally
            # Restore original resolver and package body file
            MBD._getIDCode_func[] = orig_resolver
        end
    end

    @testset "DynamicsModel constructors" begin
        # Reuse packaged body_data.xml with resolver injection
        orig_resolver = MBD._getIDCode_func[]
        MBD._getIDCode_func[] = name -> begin
            n = lowercase(strip(name))
            if n == "earth"
                return 399
            elseif n == "moon"
                return 301
            else
                return 0
            end
        end

        try
            sys = MBD.init_systemData(["Earth", "Moon"])

            # Happy path: Earth (primary), Moon (secondary)
            model = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)
            @test isa(model, MBD.CR3BPDynamicsModel)
            @test length(model.primaryData) == 2
            @test model.primaryData[1].name == "Earth"
            @test model.primaryData[2].name == "Moon"

            # Ensure the pretty-print shows the name
            s = sprint(show, model)
            @test occursin("DynamicsModel:", s) && occursin("Earth", s)

            # Empty indices -> ArgumentError
            @test_throws ArgumentError MBD.init_dynamicsModel(sys, Int64[], MBD.CR3BP)

            # Duplicate indices -> ArgumentError
            @test_throws ArgumentError MBD.init_dynamicsModel(sys, [1, 1], MBD.CR3BP)

            # Out-of-bounds index -> BoundsError
            @test_throws BoundsError MBD.init_dynamicsModel(sys, [1, 3], MBD.CR3BP)

            # Parent mismatch (secondary parent SPICEID != primary SPICEID) -> ArgumentError
            @test_throws ArgumentError MBD.init_dynamicsModel(sys, [2, 1], MBD.CR3BP)
        finally
            # Restore original resolver and package body file
            MBD._getIDCode_func[] = orig_resolver
        end
    end
end

@testset "Utilities" begin
    @testset "SPICE" begin
        @testset "getIDCode behavior" begin
            # Use the injectable function reference in the MBD SPICE utility to avoid
            # touching the external SPICE module. Save and restore the original.
            orig_getIDCode = MBD._getIDCode_func[]
            try
                # Empty name -> ArgumentError
                @test_throws ArgumentError MBD.getIDCode("")

                # Successful lookup
                MBD._getIDCode_func[] = name -> 4242
                @test MBD.getIDCode("Anything") == 4242

                # SPICE returns nothing -> ErrorException
                MBD._getIDCode_func[] = name -> nothing
                @test_throws ErrorException MBD.getIDCode("Anything")

                # SPICE throws -> ErrorException
                MBD._getIDCode_func[] = name -> throw(ErrorException("simulated failure"))
                @test_throws ErrorException MBD.getIDCode("Anything")
            finally
                MBD._getIDCode_func[] = orig_getIDCode
            end
        end
    end
end
