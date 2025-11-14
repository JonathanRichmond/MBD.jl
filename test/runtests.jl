"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 11/14/25
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
