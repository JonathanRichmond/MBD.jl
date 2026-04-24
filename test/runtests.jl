"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan LeFevre Richmond
C: 4/14/26
U: 4/24/26
"""

using MBD, Test

import LightXML, LinearAlgebra, Logging, SPICE


Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info)) # Debug/Info/Warn/Error


SPICE.furnsh("SPICEKernels/naif0012.tls", "SPICEKernels/de430.bsp")


@testset "Constructors" begin
    @testset "BodyData constructors" begin
        function clear_body_cache!()
            empty!(MBD._body_cache)
        end

        function write_xml(bodies::Vector{<:NamedTuple})::String
            path = tempname()*".xml"
            doc = LightXML.XMLDocument()
            root = LightXML.create_root(doc, "bodies")
            for b in bodies
                body_el = LightXML.new_child(root, "body")
                for (tag, val) in pairs(b)
                    child = LightXML.new_child(body_el, string(tag))
                    LightXML.set_content(child, string(val))
                end
            end
            LightXML.save_file(doc, path)
            LightXML.free(doc)

            return path
        end

        function with_mock_getIDCode(body::Function, mock_fn::Function)
            original = MBD._getIDCode_func[]
            MBD._getIDCode_func[] = mock_fn
            try
                body()
            finally
                MBD._getIDCode_func[] = original
            end
        end

        EARTH_RECORD = (
            id          = 399,
            circ_r      = 1.4959789217545033e+08,
            ecc         = 1.6735932113458880e-02,
            inc         = 2.9806674094843888e-04,
            radius      = 6.3710083666666660e+03,
            gm          = 3.9860043543609593e+05,
            raan        = 1.2204084798628088e+00,
            parentId    = 10
        )

        @testset "Input validation" begin
            mktempdir() do dir
                xml = write_xml([EARTH_RECORD])
                # Empty name throws ArgumentError
                err1 = try
                    MBD.load_bodyData("", xml)
                    nothing
                catch e
                    e
                end
                @test err1 isa ArgumentError
                @test occursin("empty", err1.msg)
                # Whitespace-only name throws ArgumentError
                for ws in ("    ", "\t", "\n", " \t\n ")
                    err2 = try
                        MBD.load_bodyData(ws, xml)
                        nothing
                    catch e
                        e
                    end
                    @test err2 isa ArgumentError
                    @test occursin("whitespace", err2.msg)
                end
                # Non-existent file throws SystemError
                err3 = try
                    MBD.load_bodyData("Earth", "/no/such/file.xml")
                    nothing
                catch e
                    e
                end
                @test err3 isa SystemError
                @test occursin("/no/such/file.xml", String(err3.prefix))
            end
        end

        @testset "Successful load" begin
            clear_body_cache!()
            mktempdir() do dir
                xml = write_xml([EARTH_RECORD])
                bd1 = MBD.load_bodyData("Earth", xml)
                # Returns BodyData
                @test bd1 isa MBD.BodyData
                # SPICE ID is correct
                @test bd1.spiceID == 399
                # Numeric fields match XML
                @test bd1.a ≈ EARTH_RECORD.circ_r
                @test bd1.e ≈ EARTH_RECORD.ecc
                @test bd1.i ≈ EARTH_RECORD.inc
                @test bd1.r ≈ EARTH_RECORD.radius
                @test bd1.μ ≈ EARTH_RECORD.gm
                @test bd1.Ω ≈ EARTH_RECORD.raan
                # Mass is derived from μ/G
                @test bd1.m ≈ EARTH_RECORD.gm/MBD.GRAVITY
                # Name is stripped input
                clear_body_cache!()
                bd2 = MBD.load_bodyData(" Earth ", xml)
                @test bd2.name == "Earth"
                # Parent SPICE ID set from XML
                @test bd1.parentSpiceID == 10
            end
        end

        @testset "Parent SPICE ID handling" begin
            clear_body_cache!()
            mktempdir() do dir
                # NaN parentId maps to 0
                record = merge(EARTH_RECORD, (parentId = "NaN",))
                xml = write_xml([record])
                bd = MBD.load_bodyData("Earth", xml)
                @test bd.parentSpiceID == MBD.UNINITIALIZED_INDEX
                clear_body_cache!()
                # nan also maps to 0
                record = merge(EARTH_RECORD, (parentId = "nan",))
                xml = write_xml([record])
                bd = MBD.load_bodyData("Earth", xml)
                @test bd.parentSpiceID == MBD.UNINITIALIZED_INDEX
                clear_body_cache!()
                # Numeric parentId is parsed correctly
                record = merge(EARTH_RECORD, (parentId = 10,))
                xml = write_xml([record])
                bd = MBD.load_bodyData("Earth", xml)
                @test bd.parentSpiceID == 10
                clear_body_cache!()
                # Unparseable parentId rethrows
                record = merge(EARTH_RECORD, (parentId = "not_a_number",))
                xml = write_xml([record])
                @test_throws Exception MBD.load_bodyData("Earth", xml)
                clear_body_cache!()
            end
        end

        @testset "XML structure errors" begin
            clear_body_cache!()
            mktempdir() do dir
                # No <body> elements throws KeyError
                path = tempname()*".xml"
                doc = LightXML.XMLDocument()
                LightXML.create_root(doc, "bodies")
                LightXML.save_file(doc, path)
                LightXML.free(doc)
                @test_throws KeyError MBD.load_bodyData("Earth", path)
                # No body matches SPICE ID throws KeyError
                mars = merge(EARTH_RECORD, (id = 499,))
                xml = write_xml([mars])
                err = try
                    MBD.load_bodyData("Earth", xml)
                    nothing
                catch e
                    e
                end
                @test err isa KeyError
                @test occursin("Earth", string(err.key))
                clear_body_cache!()
                # Missing required XML tag throws KeyError
                record = (id = 399, circ_r = 1.0, ecc = 0.0, inc = 0.0, radius = 1.0, raan = 0.0, parentId = 10)
                xml = write_xml([record])
                @test_throws KeyError MBD.load_bodyData("Earth", xml)
                clear_body_cache!()
                # Unparseable numeric field throws
                record = merge(EARTH_RECORD, (gm = "not_a_float",))
                xml = write_xml([record])
                @test_throws Exception MBD.load_bodyData("Earth", xml)
                clear_body_cache!()
                # Body with unparseable <id> is skipped; throws if no other match
                bad_id = merge(EARTH_RECORD, (id = "bad",))
                xml = write_xml([bad_id])
                @test_throws KeyError MBD.load_bodyData("Earth", xml)
                clear_body_cache!()
                # Body with unparseable <id> is skipped; valid subsequent body matches
                bad_id = merge(EARTH_RECORD, (id = "bad",))
                good_earth = EARTH_RECORD
                xml = write_xml([bad_id, good_earth])
                @test MBD.load_bodyData("Earth", xml) isa MBD.BodyData
                clear_body_cache!()
            end
        end

        @testset "SPICE lookup failure propagates" begin
            clear_body_cache!()
            mktempdir() do dir
                xml = write_xml([EARTH_RECORD])
                spice_err = ErrorException("simulated SPICE failure")
                with_mock_getIDCode((_) -> throw(spice_err)) do
                    caught = try
                        MBD.getIDCode("Earth")
                        nothing
                    catch e
                        e
                    end
                    @test caught === spice_err
                end
            end
        end

        @testset "Caching" begin
            mktempdir() do dir
                xml = write_xml([EARTH_RECORD])
                # Result stored in _body_cache after first call
                clear_body_cache!()
                MBD.load_bodyData("Earth", xml)
                @test haskey(MBD._body_cache, "Earth")
                # Stripped name is used as cache key
                clear_body_cache!()
                MBD.load_bodyData(" Earth ", xml)
                @test haskey(MBD._body_cache, "Earth")
                @test !haskey(MBD._body_cache, " Earth ")
                # Second call returns cached value without reparsing XML
                clear_body_cache!()
                MBD.load_bodyData("Earth", xml)
                fake = MBD._body_cache["Earth"]
                sentinel = MBD.BodyData(fake.a, fake.e, fake.i, fake.m, "SENTINEL", fake.parentSpiceID, fake.r, fake.spiceID, fake.μ, fake.Ω)
                MBD._body_cache["Earth"] = sentinel
                result = MBD.load_bodyData("Earth", xml)
                @test result.name == "SENTINEL"
                # Failed lookup is not cached
                clear_body_cache!()
                try
                    MBD.load_bodyData("UNKNOWN_BODY_XYZ", xml)
                catch
                end
                @test !haskey(MBD._body_cache, "UNKNOWN_BODY_XYZ")
            end
        end

        @testset "BodyData convenience constructor" begin
            clear_body_cache!()
            # BodyData(name) returns BodyData
            @test MBD.BodyData("Earth") isa MBD.BodyData
            # BodyData(name) SPICE ID matches expected
            clear_body_cache!()
            @test MBD.BodyData("Earth").spiceID == 399
            # BodyData(name) is consistent with load_bodyData on packaged file
            clear_body_cache!()
            via_constructor = MBD.BodyData("Earth")
            clear_body_cache!()
            via_loader = MBD.load_bodyData("Earth", "body_data.xml")
            @test via_constructor == via_loader
        end

        @testset "BodyData equality" begin
            clear_body_cache!()
            mktempdir() do dir
                xml = write_xml([EARTH_RECORD])
                bd1 = MBD.load_bodyData("Earth", xml)
                clear_body_cache!()
                bd2 = MBD.load_bodyData("Earth", xml)
                # Same data compares equal
                @test bd1 == bd2
                # Different SPICE ID is not equal
                mars_rec = merge(EARTH_RECORD, (id = 499,))
                xml2 = write_xml([mars_rec])
                clear_body_cache!()
                bd_mars = MBD.load_bodyData("Mars", xml2)
                @test bd1 != bd_mars
                # Different numeric field is not equal
                other_rec = merge(EARTH_RECORD, (radius = 999.0,))
                xml3 = write_xml([other_rec])
                clear_body_cache!()
                bd_other = MBD.load_bodyData("Earth", xml3)
                @test bd1 != bd_other
            end
        end

        @testset "Base.show" begin
            clear_body_cache!()
            mktempdir() do dir
                xml = write_xml([EARTH_RECORD])
                bd = MBD.load_bodyData("Earth", xml)
                # show(io, MIME, bd) does not throw
                buf = IOBuffer()
                @test_nowarn show(buf, MIME"text/plain"(), bd)
                # show output contains body name
                buf = IOBuffer()
                show(buf, MIME"text/plain"(), bd)
                @test occursin("Earth", String(take!(buf)))
                # show output contains SPICE ID
                buf = IOBuffer()
                show(buf, MIME"text/plain"(), bd)
                @test occursin("399", String(take!(buf)))
                # show(io, bd) delegates to MIME method without throwing
                buf = IOBuffer()
                @test_nowarn show(buf, bd)
            end
        end
        clear_body_cache!()
    end

    @testset "SystemData constructors" begin
        function clear_all_caches!()
            empty!(MBD._id_cache)
            empty!(MBD._body_cache)
        end

        @testset "Input validation" begin
            # Empty vector throws ArgumentError"
            err1 = try
                MBD.init_systemData(String[])
                nothing
            catch e
                e
            end
            @test err1 isa ArgumentError
            @test occursin("empty", err1.msg)
            # Vector with empty string element throws ArgumentError
            err2 = try
                MBD.init_systemData(["Earth", ""])
                nothing
            catch e
                e
            end
            @test err2 isa ArgumentError
            @test occursin("2", err2.msg)
            # Vector with whitespace-only element throws ArgumentError
            for ws in ("    ", "\t", "\n")
                err3 = try
                    MBD.init_systemData(["Earth", ws])
                    nothing
                catch e
                    e
                end
                @test err3 isa ArgumentError
                @test occursin("whitespace", err3.msg)
            end
            # Empty element at index 1 reports index 1
            err4 = try
                MBD.init_systemData([""])
                nothing
            catch e
                e
            end
            @test err4 isa ArgumentError
            @test occursin("1", err4.msg)
            # Duplicate SPICE IDs throw ArgumentError
            err5 = try
                MBD.init_systemData(["Earth", "Earth"])
                nothing
            catch e
                e
            end
            @test err5 isa ArgumentError
            @test occursin("duplicate", err5.msg)
        end

        @testset "Successful initialization" begin
            clear_all_caches!()
            # Single body returns SystemData
            @test MBD.init_systemData(["Earth"]) isa MBD.SystemData
            # Multiple bodies returns SystemData
            clear_all_caches!()
            @test MBD.init_systemData(["Earth", "Moon"]) isa MBD.SystemData
            # names field matches normalized input
            clear_all_caches!()
            sd = MBD.init_systemData([" Earth ", "Moon"])
            @test sd.names == ["Earth", "Moon"]
            # bodyData length matches names length
            @test length(sd.bodyData) == 2
            # spiceIDs length matches names length
            @test length(sd.spiceIDs) == 2
            # spiceIDs values match BodyData SPICE IDs
            @test sd.spiceIDs[1] == sd.bodyData[1].spiceID
            @test sd.spiceIDs[2] == sd.bodyData[2].spiceID
            # spiceIDs contains correct values
            @test 399 in sd.spiceIDs
            @test 301 in sd.spiceIDs
            # bodyData entries are in same order as names
            @test sd.bodyData[1].name == "Earth"
            @test sd.bodyData[2].name == "Moon"
        end

        @testset "Duplicate name handling" begin
            clear_all_caches!()
            # Duplicate names from whitespace variants still throw ArgumentError
            @test_throws ArgumentError MBD.init_systemData(["Earth", " Earth "])
        end

        @testset "BodyData loading failure propagates" begin
            clear_all_caches!()
            # Unknown body name rethrows from BodyData()
            @test_throws Exception MBD.init_systemData(["Erid"])
            # Valid body before invalid still throws on invalid
            clear_all_caches!()
            @test_throws Exception MBD.init_systemData(["Earth", "Erid"])
            # Exception is preserved (not wrapped)
            clear_all_caches!()
            err = try
                MBD.init_systemData(["Erid"])
                nothing
            catch e
                e
            end
            @test !(err isa ArgumentError)
        end

        @testset "SystemData convenience constructor" begin
            clear_all_caches!()
            # SystemData(names) returns SystemData
            @test MBD.SystemData(["Earth"]) isa MBD.SystemData
            # SystemData(names) is consistent with init_systemData
            clear_all_caches!()
            via_constructor = MBD.SystemData(["Earth", "Moon"])
            clear_all_caches!()
            via_initializer = MBD.init_systemData(["Earth", "Moon"])
            @test via_constructor == via_initializer
            # SystemData(names) propagates ArgumentError for empty input
            clear_all_caches!()
            @test_throws ArgumentError MBD.SystemData(String[])
        end

        @testset "SystemData equality" begin
            clear_all_caches!()
            sd1 = MBD.SystemData(["Earth", "Moon"])
            clear_all_caches!()
            sd2 = MBD.SystemData(["Earth", "Moon"])
            # Same data compares equal
            @test sd1 == sd2
            # Different body lists are not equal
            clear_all_caches!()
            sd_earth = MBD.SystemData(["Earth"])
            @test sd1 != sd_earth
            # Same bodies in different order are not equal
            clear_all_caches!()
            sd_order = MBD.SystemData(["Moon", "Earth"])
            @test sd1 != sd_order
        end

        @testset "Base.show" begin
            clear_all_caches!()
            sd = MBD.SystemData(["Earth", "Moon"])
            # show(io, MIME, sd) does not throw
            buf = IOBuffer()
            @test_nowarn show(buf, MIME"text/plain"(), sd)
            # show output contains all body names
            buf = IOBuffer()
            show(buf, MIME"text/plain"(), sd)
            out = String(take!(buf))
            @test occursin("Earth", out)
            @test occursin("Moon", out)
            # show output contains SPICE IDs
            buf = IOBuffer()
            show(buf, MIME"text/plain"(), sd)
            out = String(take!(buf))
            @test occursin("399", out)
            @test occursin("301", out)
            # show output contains body count
            buf = IOBuffer()
            show(buf, MIME"text/plain"(), sd)
            @test occursin("2", String(take!(buf)))
            # show singular 'body' for single-element SystemData
            clear_all_caches!()
            sd_earth = MBD.SystemData(["Earth"])
            buf = IOBuffer()
            show(buf, MIME"text/plain"(), sd_earth)
            out = String(take!(buf))
            @test occursin("1 body", out)
            @test !occursin("1 bodies", out)
            # show plural 'bodies' for multiple-element SystemData
            buf = IOBuffer()
            show(buf, MIME"text/plain"(), sd)
            @test occursin("bodies", String(take!(buf)))
            # show(io, sd) delegates to MIME method without throwing
            buf = IOBuffer()
            @test_nowarn show(buf, sd)
        end
        clear_all_caches!()
    end

#     @testset "SystemData constructors" begin
#         # Map names to the IDs present in src/body_data.xml
#         orig_resolver = MBD._getIDCode_func[]
#         MBD._getIDCode_func[] = name -> begin
#             n = lowercase(strip(name))
#             if n == "earth"
#                 return 399
#             elseif n == "moon"
#                 return 301
#             else
#                 return 0
#             end
#         end

#         try
#             sys = MBD.init_systemData(["Earth", "Moon"])
#             @test isa(sys, MBD.SystemData)
#             @test length(sys.bodyData) == 2
#             @test sys.names == ["Earth", "Moon"]
#             @test sys.SPICEIDs[1] == Int16(399)
#             @test sys.SPICEIDs[2] == Int16(301)

#             # Ensure the pretty-print shows the name
#             s = sprint(show, sys)
#             @test occursin("SystemData:", s) && occursin("Earth", s)

#             # Test equality comparison between instances
#             sys2 = MBD.init_systemData(["Earth", "Moon"])
#             @test sys == sys2
#         finally
#             # Restore original resolver and package body file
#             MBD._getIDCode_func[] = orig_resolver
#         end
#     end

#     @testset "DynamicsModel constructors" begin
#         # Reuse packaged body_data.xml with resolver injection
#         orig_resolver = MBD._getIDCode_func[]
#         MBD._getIDCode_func[] = name -> begin
#             n = lowercase(strip(name))
#             if n == "earth"
#                 return 399
#             elseif n == "moon"
#                 return 301
#             else
#                 return 0
#             end
#         end

#         try
#             sys = MBD.init_systemData(["Earth", "Moon"])

#             # Happy path: Earth (primary), Moon (secondary)
#             model = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)
#             @test isa(model, MBD.CR3BPDynamicsModel)
#             @test length(model.primaryData) == 2
#             @test model.primaryData[1].name == "Earth"
#             @test model.primaryData[2].name == "Moon"

#             # Ensure the pretty-print shows the name
#             s = sprint(show, model)
#             @test occursin("DynamicsModel:", s) && occursin("Earth", s)

#             # Test equality comparison between instances
#             model2 = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)
#             @test model == model2

#             # Empty indices -> ArgumentError
#             @test_throws ArgumentError MBD.init_dynamicsModel(sys, Int64[], MBD.CR3BP)

#             # Duplicate indices -> ArgumentError
#             @test_throws ArgumentError MBD.init_dynamicsModel(sys, [1, 1], MBD.CR3BP)

#             # Out-of-bounds index -> BoundsError
#             @test_throws BoundsError MBD.init_dynamicsModel(sys, [1, 3], MBD.CR3BP)

#             # Parent mismatch (secondary parent SPICEID != primary SPICEID) -> ArgumentError
#             @test_throws ArgumentError MBD.init_dynamicsModel(sys, [2, 1], MBD.CR3BP)
#         finally
#             # Restore original resolver and package body file
#             MBD._getIDCode_func[] = orig_resolver
#         end
#     end

#     @testset "EquationsOfMotion constructors" begin
#         # Reuse packaged body_data.xml with resolver injection
#         orig_resolver = MBD._getIDCode_func[]
#         MBD._getIDCode_func[] = name -> begin
#             n = lowercase(strip(name))
#             if n == "earth"
#                 return 399
#             elseif n == "moon"
#                 return 301
#             else
#                 return 0
#             end
#         end

#         try
#             # Build a valid CR3BP model for testing
#             sys = MBD.init_systemData(["Earth", "Moon"])
#             model = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)

#             # Happy path: Create valid CR3BPEquationsOfMotion instance
#             eom = MBD.CR3BPEquationsOfMotion(model)
#             @test isa(eom, MBD.CR3BPEquationsOfMotion)
#             @test eom.dynamicsModel === model

#             # Ensure the pretty-print shows the name
#             s = sprint(show, eom)
#             @test occursin("EquationsOfMotion:", s) && occursin("DynamicsModel", s)

#             # Test equality comparison between instances
#             eom2 = MBD.CR3BPEquationsOfMotion(model)
#             @test eom == eom2

#             # Model with wrong number of primaries (only 1 primary) -> ArgumentError
#             model_bad = MBD.CR3BPDynamicsModel([sys.bodyData[1]])
#             @test_throws ArgumentError MBD.CR3BPEquationsOfMotion(model_bad)
#         finally
#             # Restore original resolver and package body file
#             MBD._getIDCode_func[] = orig_resolver
#         end
#     end
end


# @testset "SystemData methods" begin
#     # Use resolver injection to ensure packaged body_data.xml names resolve
#     orig_resolver = MBD._getIDCode_func[]
#     MBD._getIDCode_func[] = name -> begin
#         n = lowercase(strip(name))
#         if n == "earth"
#             return 399
#         elseif n == "moon"
#             return 301
#         else
#             return 0
#         end
#     end

#     try
#         # Build a valid system
#         sys = MBD.init_systemData(["Earth", "Moon"])
#         @test isa(sys, MBD.SystemData)

#         # Test getNumPrimaries
#         @test MBD.getNumPrimaries(sys) == 2

#         # Test shallowClone
#         clone = MBD.shallowClone(sys)
#         @test isa(clone, MBD.SystemData)
#         @test clone == sys
#         @test clone.bodyData !== sys.bodyData
#         @test clone.SPICEIDs !== sys.SPICEIDs
#         @test clone.names === sys.names

#         # Error cases: inconsistent SystemData vector lengths
#         earth_bd = MBD.BodyData("Earth")
#         badSys = MBD.SystemData([earth_bd], ["Earth", "Moon"], Int16[399])
#         @test_throws ErrorException MBD.getNumPrimaries(badSys)
#         @test_throws ErrorException MBD.shallowClone(badSys)
#     finally
#         MBD._getIDCode_func[] = orig_resolver
#     end
# end


# @testset "DynamicsModel methods" begin
#     # Use resolver injection to ensure packaged body_data.xml names resolve
#     orig_resolver = MBD._getIDCode_func[]
#     MBD._getIDCode_func[] = name -> begin
#         n = lowercase(strip(name))
#         if n == "earth"
#             return 399
#         elseif n == "moon"
#             return 301
#         else
#             return 0
#         end
#     end

#     try
#         # Build a valid CR3BP model
#         sys = MBD.init_systemData(["Earth", "Moon"])
#         model = MBD.init_dynamicsModel(sys, [1, 2], MBD.CR3BP)
#         @test isa(model, MBD.CR3BPDynamicsModel)

#         # Test getNumPrimaries
#         @test MBD.getNumPrimaries(model) == 2

#         # Test appendExtraInitialConditions for CR3BP happy path - SIMPLE input
#         q0_simple = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
#         q0_simple_out = MBD.appendExtraInitialConditions(model, q0_simple, MBD.SIMPLE)
#         @test isa(q0_simple_out, Vector{Float64})
#         @test length(q0_simple_out) == 6
#         @test isapprox(q0_simple_out, q0_simple; atol=1e-12)

#         # Test appendExtraInitialConditions for CR3BP with STM
#         q0_stm = MBD.appendExtraInitialConditions(model, q0_simple, MBD.STM)
#         @test isa(q0_stm, Vector{Float64})
#         @test length(q0_stm) == 42
#         @test isapprox(q0_stm[1:6], q0_simple; atol=1e-12)
#         @test isapprox(q0_stm[8:13], zeros(Float64, 6); atol=1e-12)
#         @test isapprox(q0_stm[14], 1.0; atol=1e-12)

#         # Test appendExtraInitialConditions for CR3BP with ARCLENGTH
#         q0_arclen = MBD.appendExtraInitialConditions(model, q0_simple, MBD.ARCLENGTH)
#         @test isa(q0_arclen, Vector{Float64})
#         @test length(q0_arclen) == 43
#         @test isapprox(q0_arclen[1:6], q0_simple; atol=1e-12)
#         @test isapprox(q0_arclen[43], 0.0, atol=1e-12)

#         # Error cases: invalid input state vectors
#         @test_throws ArgumentError MBD.appendExtraInitialConditions(model, Float64[], MBD.SIMPLE)
#         @test_throws ArgumentError MBD.appendExtraInitialConditions(model, [1.0, 0.0, 0.0, 0.0, 1.0], MBD.SIMPLE)

#         # Test getCharLengths (secondary's orbital radius)
#         lstar = MBD.getCharLengths(model)
#         @test isa(lstar, Float64)
#         @test isfinite(lstar) && lstar > 0
#         @test isapprox(lstar, model.primaryData[2].a)

#         # Test getCharMasses (sum of GM/GRAVITY)
#         mstar = MBD.getCharMasses(model)
#         @test isa(mstar, Float64)
#         @test isfinite(mstar) && mstar > 0
#         expected_mstar = (model.primaryData[1].μ+model.primaryData[2].μ)/MBD.GRAVITY
#         @test isapprox(mstar, expected_mstar; rtol=1e-12)

#         # Test getCharTimes (sqrt(lstar^3/totalGM))
#         tstar = MBD.getCharTimes(model)
#         @test isa(tstar, Float64)
#         @test isfinite(tstar) && tstar > 0
#         totalGM = model.primaryData[1].μ+model.primaryData[2].μ
#         expected_tstar = sqrt(lstar^3 / totalGM)
#         @test isapprox(tstar, expected_tstar; rtol=1e-12)

#         # Test getMassRatios (μ2/(μ1+μ2))
#         μ = MBD.getMassRatios(model)
#         @test isa(μ, Float64)
#         @test 0 < μ < 1
#         expected_μ = model.primaryData[2].μ/totalGM
#         @test isapprox(μ, expected_μ; rtol=1e-12)

#         # Test getStateSize for CR3BP with all equation types
#         @test MBD.getStateSize(model, MBD.SIMPLE) == 6
#         @test MBD.getStateSize(model, MBD.STM) == 42
#         @test MBD.getStateSize(model, MBD.FULL) == 42
#         @test MBD.getStateSize(model, MBD.ARCLENGTH) == 43
#         @test MBD.getStateSize(model, MBD.MOMENTUM) == 43

#         # Verify all return Int64
#         @test isa(MBD.getStateSize(model, MBD.SIMPLE), Int64)
#         @test isa(MBD.getStateSize(model, MBD.STM), Int64)

#         # Test shallowClone
#         clone = MBD.shallowClone(model)
#         @test isa(clone, MBD.CR3BPDynamicsModel)
#         @test clone == model
#         @test clone.primaryData !== model.primaryData
#         @test clone.primaryData[1] === model.primaryData[1]

#         # Test getEquilibriumPoint for CR3BP
#         for pointID in 1:5
#             pos = MBD.getEquilibriumPoint(model, pointID)
#             @test isa(pos, Vector{Float64})
#             @test length(pos) == 3
#             @test all(isfinite, pos)
#         end

#         # Error cases: invalid pointID for getEquilibriumPoint
#         @test_throws ArgumentError MBD.getEquilibriumPoint(model, 0)
#         @test_throws ArgumentError MBD.getEquilibriumPoint(model, 6)

#         # Test getPseudopotential for CR3BP with valid position
#         pos = [1.0, 0.0, 0.0]
#         U = MBD.getPseudopotential(model, pos)
#         @test isa(U, Float64)
#         @test isfinite(U)

#         # Error cases: getPseudopotential with invalid inputs
#         @test_throws ArgumentError MBD.getPseudopotential(model, Float64[])
#         @test_throws ArgumentError MBD.getPseudopotential(model, [1.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotential(model, [NaN, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotential(model, [Inf, 0.0, 0.0])

#         # Test getPseudopotentialGradient for CR3BP with valid position
#         dU = MBD.getPseudopotentialGradient(model, pos)
#         @test isa(dU, Vector{Float64})
#         @test length(dU) == 3
#         @test all(isfinite, dU)

#         # Error cases: getPseudopotentialGradient with invalid inputs
#         @test_throws ArgumentError MBD.getPseudopotentialGradient(model, Float64[])
#         @test_throws ArgumentError MBD.getPseudopotentialGradient(model, [1.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotentialGradient(model, [NaN, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotentialGradient(model, [Inf, 0.0, 0.0])

#         # Test getPseudopotentialHessian for CR3BP with valid position
#         ddU = MBD.getPseudopotentialHessian(model, pos)
#         @test isa(ddU, Vector{Float64})
#         @test length(ddU) == 6
#         @test all(isfinite, ddU)

#         # Error cases: getPseudopotentialHessian with invalid inputs
#         @test_throws ArgumentError MBD.getPseudopotentialHessian(model, Float64[])
#         @test_throws ArgumentError MBD.getPseudopotentialHessian(model, [1.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotentialHessian(model, [NaN, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotentialHessian(model, [Inf, 0.0, 0.0])

#         # Test getEnergy for CR3BP with valid states
#         q = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]
#         JC = MBD.getEnergy(model, q)
#         @test isa(JC, Float64)
#         @test isfinite(JC)

#         # Error cases: getEnergy with invalid inputs
#         @test_throws ArgumentError MBD.getEnergy(model, Float64[])
#         @test_throws ArgumentError MBD.getEnergy(model, [1.0, 0.0, 0.0, 0.0, 1.0])
#         @test_throws ArgumentError MBD.getEnergy(model, [NaN, 0.0, 0.0, 0.0, 1.0, 0.0])
#         @test_throws ArgumentError MBD.getEnergy(model, [1.0, 0.0, 0.0, Inf, 0.1, 0.0])

#         # Test getParameterDependencies for CR3BP with FULL state vector (no parameter dependencies)
#         q_full = zeros(Float64, MBD.getStateSize(model, MBD.FULL))
#         q_full[1:6] = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]  # simple state
#         dqdparam = MBD.getParameterDependencies(model, q_full)
#         @test isa(dqdparam, Matrix{Float64})
#         @test size(dqdparam) == (6, 0)  # 6 simple states, 0 parameter sets

#         # Test with extended FULL state vector with parameter dependencies (1 parameter set)
#         q_full_extended = zeros(Float64,  MBD.getStateSize(model, MBD.FULL)+6)  # FULL + 1 param set
#         q_full_extended[1:6] = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]  # simple state
#         q_full_extended[43:48] = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]  # parameter derivatives
#         # Note: this will fail because getStateSize won't match, so we test the error instead
#         @test_throws ArgumentError MBD.getParameterDependencies(model, q_full_extended)

#         # Error cases: getParameterDependencies with invalid inputs
#         @test_throws ArgumentError MBD.getParameterDependencies(model, Float64[])
#         @test_throws ArgumentError MBD.getParameterDependencies(model, [1.0, 0.0])
#         @test_throws ArgumentError MBD.getParameterDependencies(model, [NaN; zeros(Float64, 43)])
#         @test_throws ArgumentError MBD.getParameterDependencies(model, [Inf; zeros(Float64, 43)])

#         # Test getPrimaryState for CR3BP with both primaries
#         q_primary1 = MBD.getPrimaryState(model, 1)
#         @test isa(q_primary1, Vector{Float64})
#         @test length(q_primary1) == 6
#         @test all(isfinite, q_primary1)
#         @test isapprox(q_primary1[2:6], zeros(Float64, 5); atol=1e-12)
#         expected_x1 = -μ
#         @test isapprox(q_primary1[1], expected_x1; atol=1e-12)
#         q_primary2 = MBD.getPrimaryState(model, 2)
#         @test isa(q_primary2, Vector{Float64})
#         @test length(q_primary2) == 6
#         @test all(isfinite, q_primary2)
#         @test isapprox(q_primary2[2:6], zeros(Float64, 5); atol=1e-12)
#         expected_x2 = 1 - μ
#         @test isapprox(q_primary2[1], expected_x2; atol=1e-12)

#         # Error cases: getPrimaryState with invalid inputs
#         @test_throws ArgumentError MBD.getPrimaryState(model, 0)
#         @test_throws ArgumentError MBD.getPrimaryState(model, 3)
#         @test_throws ArgumentError MBD.getPrimaryState(model, -1)

#         # Test getDistance2Primary for CR3BP with valid state
#         dist1 = MBD.getDistance2Primary(model, 1, pos)
#         @test isa(dist1, Float64)
#         @test isfinite(dist1) && dist1 >= 0
#         expected_dist1 = 1 + μ
#         @test isapprox(dist1, expected_dist1; atol=1e-12)

#         # Error cases: getDistance2Primary with invalid inputs
#         @test_throws ArgumentError MBD.getDistance2Primary(model, 1, Float64[])
#         @test_throws ArgumentError MBD.getDistance2Primary(model, 1, [1.0, 0.0])
#         @test_throws ArgumentError MBD.getDistance2Primary(model, 1, [NaN, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getDistance2Primary(model, 1, [Inf, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getDistance2Primary(model, 0, pos)
#         @test_throws ArgumentError MBD.getDistance2Primary(model, 3, pos)

#         # Test getLinearVariationState for CR3BP with 3-element variation
#         L1_pos = MBD.getEquilibriumPoint(model, 1)
#         var3 = [0.01, 0.01, 0.0]
#         q_L1, period_L1 = MBD.getLinearVariationState(model, 1, var3)
#         @test isa(q_L1, Vector{Float64})
#         @test length(q_L1) == 6
#         @test all(isfinite, q_L1)
#         @test isa(period_L1, Float64)
#         @test isfinite(period_L1) && period_L1 > 0
#         @test isapprox(q_L1[1:3], L1_pos + var3; atol=1e-10)

#         # Test getLinearVariationState with 2-element variation (z auto-appended as 0)
#         var2 = [0.01, 0.01]
#         q_L1_2d, period_L1_2d = MBD.getLinearVariationState(model, 1, var2)
#         @test isa(q_L1_2d, Vector{Float64})
#         @test length(q_L1_2d) == 6
#         @test all(isfinite, q_L1_2d)
#         @test isa(period_L1_2d, Float64)
#         @test isfinite(period_L1_2d) && period_L1_2d > 0
#         @test isapprox(q_L1_2d, q_L1; atol=1e-12)
#         @test isapprox(period_L1_2d, period_L1; atol=1e-12)

#         # Test getLinearVariationState for triangular points with Short period
#         var_tri = [0.01, 0.01, 0.0]
#         q_L4_short, period_L4_short = MBD.getLinearVariationState(model, 4, var_tri, periodType="Short")
#         @test isa(q_L4_short, Vector{Float64})
#         @test length(q_L4_short) == 6
#         @test all(isfinite, q_L4_short)
#         @test isa(period_L4_short, Float64)
#         @test isfinite(period_L4_short) && period_L4_short > 0

#         # Test getLinearVariationState for triangular points with Long period
#         q_L4_long, period_L4_long = MBD.getLinearVariationState(model, 4, var_tri, periodType="Long")
#         @test isa(q_L4_long, Vector{Float64})
#         @test length(q_L4_long) == 6
#         @test all(isfinite, q_L4_long)
#         @test isa(period_L4_long, Float64)
#         @test isfinite(period_L4_long) && period_L4_long > 0
#         @test !isapprox(period_L4_short, period_L4_long; atol=1e-2)

#         # Error cases: getLinearVariationState with invalid inputs
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 0, [0.01, 0.01])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 6, [0.01, 0.01])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 1, Float64[])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [0.01])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [0.01, 0.01, 0.01, 0.01])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [NaN, 0.01, 0.0])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 1, [0.01, Inf, 0.0])
#         @test_throws ArgumentError MBD.getLinearVariationState(model, 4, [0.01, 0.01, 0.0], periodType="Invalid")
        
#         # Test extractStateTransitionMatrix for CR3BP with valid STM state
#         q_stm = MBD.appendExtraInitialConditions(model, q, MBD.STM)
#         Φ = MBD.extractStateTransitionMatrix(model, q_stm)
#         @test isa(Φ, Matrix{Float64})
#         @test size(Φ) == (6, 6)
#         @test all(isfinite, Φ)
#         # Check that diagonal is initialized to 1 (identity for initial STM)
#         @test isapprox(LinearAlgebra.diag(Φ)[1], 1.0; atol=1e-12)

#         # Test extractStateTransitionMatrix with longer state vector (should still extract correctly)
#         q_stm_long = vcat(q_stm, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # Add extra elements
#         Φ_long = MBD.extractStateTransitionMatrix(model, q_stm_long)
#         @test isa(Φ_long, Matrix{Float64})
#         @test size(Φ_long) == (6, 6)
#         @test all(isfinite, Φ_long)
#         @test isapprox(Φ_long, Φ; atol=1e-12)

#         # Error cases: extractStateTransitionMatrix with invalid inputs
#         @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, Float64[])
#         @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, [1.0, 0.0, 0.0, 0.0, 1.0])
#         @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, [NaN; zeros(Float64, 41)])
#         @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, [Inf; zeros(Float64, 41)])
#         @test_throws ArgumentError MBD.extractStateTransitionMatrix(model, ones(Float64, 42) * Inf)

#         # Test isEpochIndependent for CR3BP
#         is_indep = MBD.isEpochIndependent(model)
#         @test isa(is_indep, Bool)
#         @test is_indep == true  # CR3BP is epoch-independent (autonomous in rotating frame)

#         # Test getEpochDependencies for CR3BP (epoch-independent)
#         q_full = zeros(Float64, MBD.getStateSize(model, MBD.FULL))
#         q_full[1:6] = [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]
#         ∂q∂E = MBD.getEpochDependencies(model, q_full)
#         @test isa(∂q∂E, Matrix{Float64})
#         @test size(∂q∂E) == (6, 0)  # Empty matrix for epoch-independent CR3BP
#         @test all(isfinite, ∂q∂E)

#         # Error cases: getEpochDependencies with invalid inputs
#         @test_throws ArgumentError MBD.getEpochDependencies(model, Float64[])
#         @test_throws ArgumentError MBD.getEpochDependencies(model, [1.0, 0.0, 0.0, 0.0, 1.0])
#         @test_throws ArgumentError MBD.getEpochDependencies(model, [NaN; zeros(Float64, 41)])
#         @test_throws ArgumentError MBD.getEpochDependencies(model, [Inf; zeros(Float64, 41)])
#         @test_throws ArgumentError MBD.getEpochDependencies(model, ones(Float64, 45))

#         # Test getEquationsOfMotion for CR3BP with valid model
#         eom = MBD.getEquationsOfMotion(model)
#         @test isa(eom, MBD.CR3BPEquationsOfMotion)
#         @test eom.dynamicsModel === model
#         @test isa(eom.dynamicsModel, MBD.CR3BPDynamicsModel)
        
#         # Error cases: abstract methods throw on non-CR3BP model
#         struct testModel <: MBD.AbstractDynamicsModel end
#         tm = testModel()
#         @test_throws ErrorException MBD.appendExtraInitialConditions(tm, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0], MBD.SIMPLE)
#         @test_throws ErrorException MBD.extractStateTransitionMatrix(tm, ones(Float64, 42))
#         @test_throws ErrorException MBD.isEpochIndependent(tm)
#         @test_throws ErrorException MBD.getCharLengths(tm)
#         @test_throws ErrorException MBD.getCharMasses(tm)
#         @test_throws ErrorException MBD.getCharTimes(tm)
#         @test_throws ErrorException MBD.getDistance2Primary(tm, 1, [1.0, 0.0, 0.0])
#         @test_throws ErrorException MBD.getEquationsOfMotion(tm)
#         @test_throws ErrorException MBD.getEquilibriumPoint(tm, 1)
#         @test_throws ErrorException MBD.getMassRatios(tm)
#         @test_throws ErrorException MBD.getPrimaryState(tm, 1)
#         @test_throws ErrorException MBD.getPseudopotential(tm, [1.0, 0.0, 0.0])
#         @test_throws ErrorException MBD.getPseudopotentialGradient(tm, [1.0, 0.0, 0.0])
#         @test_throws ErrorException MBD.getPseudopotentialHessian(tm, [1.0, 0.0, 0.0])
#         @test_throws ErrorException MBD.getStateSize(tm, MBD.SIMPLE)
#         @test_throws ErrorException MBD.shallowClone(tm)

#         # Error cases: getNumPrimaries with no primaryData field
#         @test_throws ErrorException MBD.getNumPrimaries(tm)

#         # Error cases: CR3BP methods with non-CR3BP model
#         @test_throws MethodError MBD.getLinearVariationState(tm, 1, [0.01, 0.01, 0.0])

#         # Error cases: CR3BP methods with wrong primary count
#         # Create a malformed model with only 1 body (for testing purposes)
#         bad_model = MBD.CR3BPDynamicsModel([model.primaryData[1]])
#         @test_throws ArgumentError MBD.appendExtraInitialConditions(bad_model, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0], MBD.SIMPLE)
#         @test_throws ArgumentError MBD.extractStateTransitionMatrix(bad_model, ones(Float64, 42))
#         @test_throws ArgumentError MBD.isEpochIndependent(bad_model)
#         @test_throws ArgumentError MBD.getCharLengths(bad_model)
#         @test_throws ArgumentError MBD.getCharMasses(bad_model)
#         @test_throws ArgumentError MBD.getCharTimes(bad_model)
#         @test_throws ArgumentError MBD.getDistance2Primary(bad_model, 1, [1.0, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getEnergy(bad_model, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
#         @test_throws ArgumentError MBD.getEquationsOfMotion(bad_model)
#         @test_throws ArgumentError MBD.getEquilibriumPoint(bad_model, 1)
#         @test_throws ArgumentError MBD.getLinearVariationState(bad_model, 1, [0.01, 0.01, 0.0])
#         @test_throws ArgumentError MBD.getMassRatios(bad_model)
#         @test_throws ArgumentError MBD.getPrimaryState(bad_model, 1)
#         @test_throws ArgumentError MBD.getPseudopotential(bad_model, [1.0, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotentialGradient(bad_model, [1.0, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getPseudopotentialHessian(bad_model, [1.0, 0.0, 0.0])
#         @test_throws ArgumentError MBD.getStateSize(bad_model, MBD.SIMPLE)
#     finally
#         MBD._getIDCode_func[] = orig_resolver
#     end
# end

@testset "Utilities" begin
    @testset "SPICE" begin
        @testset "getIDCode behavior" begin
            function with_mock_getIDCode(body::Function, mock_fn::Function)
                original = MBD._getIDCode_func[]
                MBD._getIDCode_func[] = mock_fn
                try
                    body()
                finally
                    MBD._getIDCode_func[] = original
                end
            end

            function clear_id_cache!()
                empty!(MBD._id_cache)
            end

            @testset "Input validation" begin
                # Empty string throws ArgumentError
                err1 = try
                    MBD.getIDCode("")
                    nothing
                catch e
                    e
                end
                @test err1 isa ArgumentError
                @test occursin("empty", err1.msg)
                # Whitespace-only string throws ArgumentError
                for ws in ("    ", "\t", "\n", " \t\n ")
                    err2 = try
                        MBD.getIDCode(ws)
                        nothing
                    catch e
                        e
                    end
                    @test err2 isa ArgumentError
                    @test occursin("whitespace", err2.msg)
                end
            end

            @testset "Successful SPICE resolution" begin
                clear_id_cache!()
                # Known body returns correct ID
                code = MBD.getIDCode("Earth")
                @test code == 399
                @test code isa Int64
                # Case/whitespace is stripped before lookup
                clear_id_cache!()
                @test MBD.getIDCode(" EARTH ") == 399
                clear_id_cache!()
            end

            @testset "Caching" begin
                clear_id_cache!()
                # Result is stored in _id_cache after first call
                code = MBD.getIDCode("Earth")
                @test haskey(MBD._id_cache, "Earth")
                @test MBD._id_cache["Earth"] == 399
                # Second call return cached value without hitting SPICE
                clear_id_cache!()
                callCount = Ref(0)
                mock = function(name::String)
                    normalized = strip(name)
                    if haskey(MBD._id_cache, normalized)
                        return MBD._id_cache[normalized]
                    end
                    local code::Int
                    callCount[] += 1
                    code = SPICE.bods2c(normalized)
                    MBD._id_cache[normalized] = code
                    return code
                end
                with_mock_getIDCode(mock) do 
                    code1 = getIDCode("Earth")
                    spiceCallsBefore = callCount[]
                    code2 = getIDCode("Earth")
                    @test callCount[] == spiceCallsBefore
                end
                # Whitespace-normalized name is cached under stripped key
                clear_id_cache!()
                code = MBD.getIDCode(" Earth ")
                @test haskey(MBD._id_cache, "Earth")
                @test !haskey(MBD._id_cache, " Earth ")
                # Cache hit returns same value as original resolution
                clear_id_cache!()
                code1 = MBD.getIDCode("Earth")
                code2 = MBD.getIDCode("Earth")
                @test code1 == code2
                clear_id_cache!()
            end

            @testset "Unknown body name" begin
                clear_id_cache!()
                # Unrecognized body throws KeyError
                err = try
                    MBD.getIDCode("Erid")
                    nothing
                catch e
                    e
                end
                @test err isa KeyError
                @test occursin("Erid", String(err.key))
                # Unrecognized body is not stored in cache
                @test !haskey(MBD._id_cache, "Erid")
                clear_id_cache!()
            end

            @testset "SPICE exception propagation" begin
                clear_id_cache!()
                # Arbitrary SPICE error is re-thrown
                mock = function(name::String)
                    normalized = strip(name)
                    if haskey(MBD._id_cache, normalized)
                        return MBD._id_cache[normalized]
                    end
                    if normalized == "Mars"
                        throw(DomainError("simulated SPICE kernel failure"))
                    end
                end
                with_mock_getIDCode(mock) do
                    caught = try
                        MBD.getIDCode("Mars")
                        nothing
                    catch e
                        e
                    end
                    @test caught === DomainError("simulated SPICE kernel failure")
                    # Exception type is preserved
                    @test_throws DomainError MBD.getIDCode("Mars")
                end
            end
        end
    end
end


SPICE.kclear()
