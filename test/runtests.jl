"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan LeFevre Richmond
C: 12/12/25
U: 1/15/26
"""

using MBD, Test

import LightXML, LinearAlgebra, Logging, SPICE


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

            # Test equality comparison between instances
            bd2 = MBD.load_bodyData("TestBody", tmpfile)
            @test bd == bd2
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

            # Test equality comparison between instances
            sys2 = MBD.init_systemData(["Earth", "Moon"])
            @test sys == sys2
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

            # Test equality comparison between instances
            model2 = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)
            @test model == model2

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

    @testset "EquationsOfMotion constructors" begin
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
            # Build a valid CR3BP model for testing
            sys = MBD.init_systemData(["Earth", "Moon"])
            model = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)

            # Happy path: Create valid CR3BPEquationsOfMotion instance
            eom = MBD.CR3BPEquationsOfMotion(model)
            @test isa(eom, MBD.CR3BPEquationsOfMotion)
            @test eom.dynamicsModel === model

            # Ensure the pretty-print shows the name
            s = sprint(show, eom)
            @test occursin("EquationsOfMotion:", s) && occursin("DynamicsModel", s)

            # Test equality comparison between instances
            eom2 = MBD.CR3BPEquationsOfMotion(model)
            @test eom == eom2

            # Model with wrong number of primaries (only 1 primary) -> ArgumentError
            model_bad = MBD.CR3BPDynamicsModel([sys.bodyData[1]])
            @test_throws ArgumentError MBD.CR3BPEquationsOfMotion(model_bad)
        finally
            # Restore original resolver and package body file
            MBD._getIDCode_func[] = orig_resolver
        end
    end
end


@testset "SystemData methods" begin
    # Use resolver injection to ensure packaged body_data.xml names resolve
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
        # Build a valid system
        sys = MBD.init_systemData(["Earth", "Moon"])
        @test isa(sys, MBD.SystemData)

        # Test getNumPrimaries
        @test MBD.getNumPrimaries(sys) == 2

        # Test shallowClone
        clone = MBD.shallowClone(sys)
        @test isa(clone, MBD.SystemData)
        @test clone == sys
        @test clone.bodyData !== sys.bodyData
        @test clone.SPICEIDs !== sys.SPICEIDs
        @test clone.names === sys.names

        # Error cases: inconsistent SystemData vector lengths
        earth_bd = MBD.BodyData("Earth")
        badSys = MBD.SystemData([earth_bd], ["Earth", "Moon"], Int16[399])
        @test_throws ErrorException MBD.getNumPrimaries(badSys)
        @test_throws ErrorException MBD.shallowClone(badSys)
    finally
        MBD._getIDCode_func[] = orig_resolver
    end
end


@testset "DynamicsModel methods" begin
    # Use resolver injection to ensure packaged body_data.xml names resolve
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
        # Build a valid CR3BP model
        sys = MBD.init_systemData(["Earth", "Moon"])
        model = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)
        @test isa(model, MBD.CR3BPDynamicsModel)

        # Test getNumPrimaries
        @test MBD.getNumPrimaries(model) == 2

        # Test appendExtraInitialConditions for CR3BP happy path - SIMPLE input
        q0_simple = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        q0_simple_out = MBD.appendExtraInitialConditions(model, q0_simple, MBD.SIMPLE)
        @test isa(q0_simple_out, Vector{Float64})
        @test length(q0_simple_out) == 6
        @test isapprox(q0_simple_out, q0_simple; atol=1e-12)

        # Test appendExtraInitialConditions for CR3BP with STM
        q0_stm = MBD.appendExtraInitialConditions(model, q0_simple, MBD.STM)
        @test isa(q0_stm, Vector{Float64})
        @test length(q0_stm) == 42
        @test isapprox(q0_stm[1:6], q0_simple; atol=1e-12)
        @test isapprox(q0_stm[8:13], zeros(Float64, 6); atol=1e-12)
        @test isapprox(q0_stm[14], 1.0; atol=1e-12)

        # Test appendExtraInitialConditions for CR3BP with ARCLENGTH
        q0_arclen = MBD.appendExtraInitialConditions(model, q0_simple, MBD.ARCLENGTH)
        @test isa(q0_arclen, Vector{Float64})
        @test length(q0_arclen) == 43
        @test isapprox(q0_arclen[1:6], q0_simple; atol=1e-12)
        @test isapprox(q0_arclen[43], 0.0, atol=1e-12)

        # Error cases: invalid input state vectors
        @test_throws ArgumentError MBD.appendExtraInitialConditions(model, Float64[], MBD.SIMPLE)
        @test_throws ArgumentError MBD.appendExtraInitialConditions(model, [1.0, 0.0, 0.0, 0.0, 1.0], MBD.SIMPLE)

        # Test getCharLengths (secondary's orbital radius)
        lstar = MBD.getCharLengths(model)
        @test isa(lstar, Float64)
        @test isfinite(lstar) && lstar > 0
        @test isapprox(lstar, model.primaryData[2].a)

        # Test getCharMasses (sum of GM/GRAVITY)
        mstar = MBD.getCharMasses(model)
        @test isa(mstar, Float64)
        @test isfinite(mstar) && mstar > 0
        expected_mstar = (model.primaryData[1].μ+model.primaryData[2].μ)/MBD.GRAVITY
        @test isapprox(mstar, expected_mstar; rtol=1e-12)

        # Test getCharTimes (sqrt(lstar^3/totalGM))
        tstar = MBD.getCharTimes(model)
        @test isa(tstar, Float64)
        @test isfinite(tstar) && tstar > 0
        totalGM = model.primaryData[1].μ+model.primaryData[2].μ
        expected_tstar = sqrt(lstar^3 / totalGM)
        @test isapprox(tstar, expected_tstar; rtol=1e-12)

        # Test getMassRatios (μ2/(μ1+μ2))
        μ = MBD.getMassRatios(model)
        @test isa(μ, Float64)
        @test 0 < μ < 1
        expected_μ = model.primaryData[2].μ/totalGM
        @test isapprox(μ, expected_μ; rtol=1e-12)

        # Test getStateSize for CR3BP with all equation types
        @test MBD.getStateSize(model, MBD.SIMPLE) == 6
        @test MBD.getStateSize(model, MBD.STM) == 42
        @test MBD.getStateSize(model, MBD.FULL) == 42
        @test MBD.getStateSize(model, MBD.ARCLENGTH) == 43
        @test MBD.getStateSize(model, MBD.MOMENTUM) == 43

        # Verify all return Int64
        @test isa(MBD.getStateSize(model, MBD.SIMPLE), Int64)
        @test isa(MBD.getStateSize(model, MBD.STM), Int64)

        # Test shallowClone
        clone = MBD.shallowClone(model)
        @test isa(clone, MBD.CR3BPDynamicsModel)
        @test clone == model
        @test clone.primaryData !== model.primaryData
        @test clone.primaryData[1] === model.primaryData[1]

        # Test getEquilibriumPoint for CR3BP
        for pointID in 1:5
            pos = MBD.getEquilibriumPoint(model, pointID)
            @test isa(pos, Vector{Float64})
            @test length(pos) == 3
            @test all(isfinite, pos)
        end

        # Error cases: invalid pointID for getEquilibriumPoint
        @test_throws ArgumentError MBD.getEquilibriumPoint(model, 0)
        @test_throws ArgumentError MBD.getEquilibriumPoint(model, 6)

        # Test getPseudopotential for CR3BP with valid position
        pos = [1.0, 0.0, 0.0]
        U = MBD.getPseudopotential(model, pos)
        @test isa(U, Float64)
        @test isfinite(U)

        # Error cases: getPseudopotential with invalid inputs
        @test_throws ArgumentError MBD.getPseudopotential(model, Float64[])
        @test_throws ArgumentError MBD.getPseudopotential(model, [1.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotential(model, [NaN, 0.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotential(model, [Inf, 0.0, 0.0])

        # Test getPseudopotentialGradient for CR3BP with valid position
        dU = MBD.getPseudopotentialGradient(model, pos)
        @test isa(dU, Vector{Float64})
        @test length(dU) == 3
        @test all(isfinite, dU)

        # Error cases: getPseudopotentialGradient with invalid inputs
        @test_throws ArgumentError MBD.getPseudopotentialGradient(model, Float64[])
        @test_throws ArgumentError MBD.getPseudopotentialGradient(model, [1.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotentialGradient(model, [NaN, 0.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotentialGradient(model, [Inf, 0.0, 0.0])

        # Test getPseudopotentialHessian for CR3BP with valid position
        ddU = MBD.getPseudopotentialHessian(model, pos)
        @test isa(ddU, Vector{Float64})
        @test length(ddU) == 6
        @test all(isfinite, ddU)

        # Error cases: getPseudopotentialHessian with invalid inputs
        @test_throws ArgumentError MBD.getPseudopotentialHessian(model, Float64[])
        @test_throws ArgumentError MBD.getPseudopotentialHessian(model, [1.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotentialHessian(model, [NaN, 0.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotentialHessian(model, [Inf, 0.0, 0.0])

        # Test getEnergy for CR3BP with valid states
        q = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]
        JC = MBD.getEnergy(model, q)
        @test isa(JC, Float64)
        @test isfinite(JC)

        # Error cases: getEnergy with invalid inputs
        @test_throws ArgumentError MBD.getEnergy(model, Float64[])
        @test_throws ArgumentError MBD.getEnergy(model, [1.0, 0.0, 0.0, 0.0, 1.0])
        @test_throws ArgumentError MBD.getEnergy(model, [NaN, 0.0, 0.0, 0.0, 1.0, 0.0])
        @test_throws ArgumentError MBD.getEnergy(model, [1.0, 0.0, 0.0, Inf, 0.1, 0.0])

        # Test getParameterDependencies for CR3BP with FULL state vector (no parameter dependencies)
        q_full = zeros(Float64, MBD.getStateSize(model, MBD.FULL))
        q_full[1:6] = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]  # simple state
        dqdparam = MBD.getParameterDependencies(model, q_full)
        @test isa(dqdparam, Matrix{Float64})
        @test size(dqdparam) == (6, 0)  # 6 simple states, 0 parameter sets

        # Test with extended FULL state vector with parameter dependencies (1 parameter set)
        q_full_extended = zeros(Float64,  MBD.getStateSize(model, MBD.FULL)+6)  # FULL + 1 param set
        q_full_extended[1:6] = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]  # simple state
        q_full_extended[43:48] = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]  # parameter derivatives
        # Note: this will fail because getStateSize won't match, so we test the error instead
        @test_throws ArgumentError MBD.getParameterDependencies(model, q_full_extended)

        # Error cases: getParameterDependencies with invalid inputs
        @test_throws ArgumentError MBD.getParameterDependencies(model, Float64[])
        @test_throws ArgumentError MBD.getParameterDependencies(model, [1.0, 0.0])
        @test_throws ArgumentError MBD.getParameterDependencies(model, [NaN; zeros(Float64, 43)])
        @test_throws ArgumentError MBD.getParameterDependencies(model, [Inf; zeros(Float64, 43)])

        # Test getPrimaryState for CR3BP with both primaries
        q_primary1 = MBD.getPrimaryState(model, 1)
        @test isa(q_primary1, Vector{Float64})
        @test length(q_primary1) == 6
        @test all(isfinite, q_primary1)
        @test isapprox(q_primary1[2:6], zeros(Float64, 5); atol=1e-12)
        expected_x1 = -μ
        @test isapprox(q_primary1[1], expected_x1; atol=1e-12)
        q_primary2 = MBD.getPrimaryState(model, 2)
        @test isa(q_primary2, Vector{Float64})
        @test length(q_primary2) == 6
        @test all(isfinite, q_primary2)
        @test isapprox(q_primary2[2:6], zeros(Float64, 5); atol=1e-12)
        expected_x2 = 1 - μ
        @test isapprox(q_primary2[1], expected_x2; atol=1e-12)

        # Error cases: getPrimaryState with invalid inputs
        @test_throws ArgumentError MBD.getPrimaryState(model, 0)
        @test_throws ArgumentError MBD.getPrimaryState(model, 3)
        @test_throws ArgumentError MBD.getPrimaryState(model, -1)

        # Test getDistance2Primary for CR3BP with valid state
        dist1 = MBD.getDistance2Primary(model, 1, pos)
        @test isa(dist1, Float64)
        @test isfinite(dist1) && dist1 >= 0
        expected_dist1 = 1 + μ
        @test isapprox(dist1, expected_dist1; atol=1e-12)

        # Error cases: getDistance2Primary with invalid inputs
        @test_throws ArgumentError MBD.getDistance2Primary(model, 1, Float64[])
        @test_throws ArgumentError MBD.getDistance2Primary(model, 1, [1.0, 0.0])
        @test_throws ArgumentError MBD.getDistance2Primary(model, 1, [NaN, 0.0, 0.0])
        @test_throws ArgumentError MBD.getDistance2Primary(model, 1, [Inf, 0.0, 0.0])
        @test_throws ArgumentError MBD.getDistance2Primary(model, 0, pos)
        @test_throws ArgumentError MBD.getDistance2Primary(model, 3, pos)

        # Test getLinearVariationState for CR3BP with 3-element variation
        L1_pos = MBD.getEquilibriumPoint(model, 1)
        var3 = [0.01, 0.01, 0.0]
        q_L1, period_L1 = MBD.getLinearVariationState(model, 1, var3)
        @test isa(q_L1, Vector{Float64})
        @test length(q_L1) == 6
        @test all(isfinite, q_L1)
        @test isa(period_L1, Float64)
        @test isfinite(period_L1) && period_L1 > 0
        @test isapprox(q_L1[1:3], L1_pos + var3; atol=1e-10)

        # Test getLinearVariationState with 2-element variation (z auto-appended as 0)
        var2 = [0.01, 0.01]
        q_L1_2d, period_L1_2d = MBD.getLinearVariationState(model, 1, var2)
        @test isa(q_L1_2d, Vector{Float64})
        @test length(q_L1_2d) == 6
        @test all(isfinite, q_L1_2d)
        @test isa(period_L1_2d, Float64)
        @test isfinite(period_L1_2d) && period_L1_2d > 0
        @test isapprox(q_L1_2d, q_L1; atol=1e-12)
        @test isapprox(period_L1_2d, period_L1; atol=1e-12)

        # Test getLinearVariationState for triangular points with Short period
        var_tri = [0.01, 0.01, 0.0]
        q_L4_short, period_L4_short = MBD.getLinearVariationState(model, 4, var_tri, periodType="Short")
        @test isa(q_L4_short, Vector{Float64})
        @test length(q_L4_short) == 6
        @test all(isfinite, q_L4_short)
        @test isa(period_L4_short, Float64)
        @test isfinite(period_L4_short) && period_L4_short > 0

        # Test getLinearVariationState for triangular points with Long period
        q_L4_long, period_L4_long = MBD.getLinearVariationState(model, 4, var_tri, periodType="Long")
        @test isa(q_L4_long, Vector{Float64})
        @test length(q_L4_long) == 6
        @test all(isfinite, q_L4_long)
        @test isa(period_L4_long, Float64)
        @test isfinite(period_L4_long) && period_L4_long > 0
        @test !isapprox(period_L4_short, period_L4_long; atol=1e-2)

        # Error cases: getLinearVariationState with invalid inputs
        @test_throws ArgumentError MBD.getLinearVariationState(model, 0, [0.01, 0.01])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 6, [0.01, 0.01])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 1, Float64[])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [0.01])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [0.01, 0.01, 0.01, 0.01])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [NaN, 0.01, 0.0])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [0.01, Inf, 0.0])
        @test_throws ArgumentError MBD.getLinearVariationState(model, 4, [0.01, 0.01, 0.0], periodType="Invalid")
        
        # Test extractStateTransitionMatrix for CR3BP with valid STM state
        q_stm = MBD.appendExtraInitialConditions(model, q, MBD.STM)
        Φ = MBD.extractStateTransitionMatrix(model, q_stm)
        @test isa(Φ, Matrix{Float64})
        @test size(Φ) == (6, 6)
        @test all(isfinite, Φ)
        # Check that diagonal is initialized to 1 (identity for initial STM)
        @test isapprox(LinearAlgebra.diag(Φ)[1], 1.0; atol=1e-12)

        # Test extractStateTransitionMatrix with longer state vector (should still extract correctly)
        q_stm_long = vcat(q_stm, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # Add extra elements
        Φ_long = MBD.extractStateTransitionMatrix(model, q_stm_long)
        @test isa(Φ_long, Matrix{Float64})
        @test size(Φ_long) == (6, 6)
        @test all(isfinite, Φ_long)
        @test isapprox(Φ_long, Φ; atol=1e-12)

        # Error cases: extractStateTransitionMatrix with invalid inputs
        @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, Float64[])
        @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, [1.0, 0.0, 0.0, 0.0, 1.0])
        @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, [NaN; zeros(Float64, 41)])
        @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, [Inf; zeros(Float64, 41)])
        @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, ones(Float64, 42) * Inf)

        # Test isEpochIndependent for CR3BP
        is_indep = MBD.isEpochIndependent(model)
        @test isa(is_indep, Bool)
        @test is_indep == true  # CR3BP is epoch-independent (autonomous in rotating frame)

        # Test getEpochDependencies for CR3BP (epoch-independent)
        q_full = zeros(Float64, MBD.getStateSize(model, MBD.FULL))
        q_full[1:6] = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]
        ∂q∂E = MBD.getEpochDependencies(model, q_full)
        @test isa(∂q∂E, Matrix{Float64})
        @test size(∂q∂E) == (6, 0)  # Empty matrix for epoch-independent CR3BP
        @test all(isfinite, ∂q∂E)

        # Error cases: getEpochDependencies with invalid inputs
        @test_throws ArgumentError MBD.getEpochDependencies(model, Float64[])
        @test_throws ArgumentError MBD.getEpochDependencies(model, [1.0, 0.0, 0.0, 0.0, 1.0])
        @test_throws ArgumentError MBD.getEpochDependencies(model, [NaN; zeros(Float64, 41)])
        @test_throws ArgumentError MBD.getEpochDependencies(model, [Inf; zeros(Float64, 41)])
        @test_throws ArgumentError MBD.getEpochDependencies(model, ones(Float64, 45))
        
        # Error cases: abstract methods throw on non-CR3BP model
        struct testModel <: MBD.AbstractDynamicsModel end
        tm = testModel()
        @test_throws ErrorException MBD.appendExtraInitialConditions(tm, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0], MBD.SIMPLE)
        @test_throws ErrorException MBD.extractStateTransitionMatrix(tm, ones(Float64, 42))
        @test_throws ErrorException MBD.isEpochIndependent(tm)
        @test_throws ErrorException MBD.getCharLengths(tm)
        @test_throws ErrorException MBD.getCharMasses(tm)
        @test_throws ErrorException MBD.getCharTimes(tm)
        @test_throws ErrorException MBD.getDistance2Primary(tm, 1, [1.0, 0.0, 0.0])
        @test_throws ErrorException MBD.getEquilibriumPoint(tm, 1)
        @test_throws ErrorException MBD.getMassRatios(tm)
        @test_throws ErrorException MBD.getPrimaryState(tm, 1)
        @test_throws ErrorException MBD.getPseudopotential(tm, [1.0, 0.0, 0.0])
        @test_throws ErrorException MBD.getPseudopotentialGradient(tm, [1.0, 0.0, 0.0])
        @test_throws ErrorException MBD.getPseudopotentialHessian(tm, [1.0, 0.0, 0.0])
        @test_throws ErrorException MBD.getStateSize(tm, MBD.SIMPLE)
        @test_throws ErrorException MBD.shallowClone(tm)

        # Error cases: getNumPrimaries with no primaryData field
        @test_throws ErrorException MBD.getNumPrimaries(tm)

        # Error cases: CR3BP methods with non-CR3BP model
        @test_throws MethodError MBD.getLinearVariationState(tm, 1, [0.01, 0.01, 0.0])

        # Error cases: CR3BP methods with wrong primary count
        # Create a malformed model with only 1 body (for testing purposes)
        bad_model = MBD.CR3BPDynamicsModel([model.primaryData[1]])
        @test_throws ArgumentError MBD.appendExtraInitialConditions(bad_model, q0_simple, MBD.SIMPLE)
        @test_throws ArgumentError MBD.extractStateTransitionMatrix(bad_model, ones(Float64, 42))
        @test_throws ArgumentError MBD.isEpochIndependent(bad_model)
        @test_throws ArgumentError MBD.getCharLengths(bad_model)
        @test_throws ArgumentError MBD.getCharMasses(bad_model)
        @test_throws ArgumentError MBD.getCharTimes(bad_model)
        @test_throws ArgumentError MBD.getDistance2Primary(bad_model, 1, [1.0, 0.0, 0.0])
        @test_throws ArgumentError MBD.getEquilibriumPoint(bad_model, 1)
        @test_throws ArgumentError MBD.getLinearVariationState(bad_model, 1, [0.01, 0.01, 0.0])
        @test_throws ArgumentError MBD.getMassRatios(bad_model)
        @test_throws ArgumentError MBD.getPrimaryState(bad_model, 1)
        @test_throws ArgumentError MBD.getPseudopotential(bad_model, [1.0, 0.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotentialGradient(bad_model, [1.0, 0.0, 0.0])
        @test_throws ArgumentError MBD.getPseudopotentialHessian(bad_model, [1.0, 0.0, 0.0])
        @test_throws ArgumentError MBD.getStateSize(bad_model, MBD.SIMPLE)
    finally
        MBD._getIDCode_func[] = orig_resolver
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
