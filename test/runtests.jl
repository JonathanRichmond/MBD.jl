"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan Richmond
C: 11/14/25
"""

using Logging, MBD, Test

global_logger(ConsoleLogger(stderr, Logging.Info)) # Debug/Info/Warn/Error

@testset "Files" begin
    @test isfile("../src/body_data.xml")
end
