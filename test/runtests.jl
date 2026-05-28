"""
Multi-Body Dynamics astrodynamics package tests

Author: Jonathan LeFevre Richmond
C: 4/14/26
U: 5/28/26
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
                @test isequal(bd1, bd2)
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
                out = String(take!(buf))
                @test occursin("Earth", out)
                # show output contains SPICE ID
                @test occursin("399", out)
                # show(io, bd) delegates to MIME method without throwing
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
            @test isequal(sd1, sd2)
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
            out = String(take!(buf))
            @test occursin("Earth", out)
            @test occursin("Moon", out)
            # show output contains SPICE IDs
            @test occursin("399", out)
            @test occursin("301", out)
            # show output contains body count
            @test occursin("2", out)
            # show singular 'body' for single-element SystemData
            clear_all_caches!()
            sd_earth = MBD.SystemData(["Earth"])
            show(buf, MIME"text/plain"(), sd_earth)
            out = String(take!(buf))
            @test occursin("1 body", out)
            @test !occursin("1 bodies", out)
            # show plural 'bodies' for multiple-element SystemData
            show(buf, MIME"text/plain"(), sd)
            @test occursin("bodies", String(take!(buf)))
            # show(io, sd) delegates to MIME method without throwing
            @test_nowarn show(buf, sd)
        end
        clear_all_caches!()
    end

    @testset "DynamicsModel constructors" begin
        function clear_all_caches!()
            empty!(MBD._id_cache)
            empty!(MBD._body_cache)
        end

        EARTH_MOON_SYSTEM = let
            clear_all_caches!()
            MBD.SystemData(["Earth", "Moon"])
        end

        SUN_EARTH_MOON_SYSTEM = let
            clear_all_caches!()
            MBD.SystemData(["Sun", "Earth", "Moon"])
        end

        EARTH_MOON_SUN_SYSTEM = let
            clear_all_caches!()
            MBD.SystemData(["Earth", "Moon", "Sun"])
        end

        @testset "Input validation" begin
            # Empty indices throws ArgumentError
            err1 = try
                MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, Vector{Int64}(), MBD.CR3BPDynamicsModel)
                nothing
            catch e
                e
            end
            @test err1 isa ArgumentError
            @test occursin("empty", err1.msg)
            # Index of 0 throws BoundsError
            @test_throws BoundsError MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [0, 1], MBD.CR3BPDynamicsModel)
            # Index beyond length throws BoundsError
            n = length(EARTH_MOON_SYSTEM.bodyData)
            @test_throws BoundsError MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [1, n+1], MBD.CR3BPDynamicsModel)
            # Negative index throws BoundsError
            @test_throws BoundsError MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [-1, 1], MBD.CR3BPDynamicsModel)
            # Duplicate indices throw ArgumentError
            err2 = try
                MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [1, 1], MBD.CR3BPDynamicsModel)
                nothing
            catch e
                e
            end
            @test err2 isa ArgumentError
            @test occursin("duplicate", err2.msg)
            # BoundsError is checked before duplicate check
            @test_throws BoundsError MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [n+1, n+1], MBD.CR3BPDynamicsModel)
        end

        @testset "CR3BP validation" begin
            # Fewer than 2 bodies throws ArgumentError
            err1 = try
                MBD.build_dynamicsModel(MBD.CR3BPDynamicsModel, [EARTH_MOON_SYSTEM.bodyData[1]])
                nothing
            catch e
                e
            end
            @test err1 isa ArgumentError
            @test occursin("2", err1.msg)
            # More than 2 bodies throws ArgumentError
            @test_throws ArgumentError MBD.build_dynamicsModel(MBD.CR3BPDynamicsModel, SUN_EARTH_MOON_SYSTEM.bodyData)
            # Reversed body order (secondary not child of primary) throws ArgumentError
            err2 = try
                MBD.build_dynamicsModel(MBD.CR3BPDynamicsModel, reverse(EARTH_MOON_SYSTEM.bodyData))
                nothing
            catch e
                e
            end
            @test err2 isa ArgumentError
            @test occursin("parent", err2.msg)
            # Unrelated bodies throws ArgumentError
            clear_all_caches!()
            sd = MBD.SystemData(["Earth", "Mars"])
            @test_throws ArgumentError MBD.build_dynamicsModel(MBD.CR3BPDynamicsModel, sd.bodyData)
        end

        @testset "Successful initialization" begin
            # init_dynamicsModel returns CR3BPDynamicsModel
            dm = MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [1, 2], MBD.CR3BPDynamicsModel)
            @test dm isa MBD.CR3BPDynamicsModel
            # primaryData has exactly 2 entries
            @test length(dm.primaryData) == 2
            # primaryData entries match selected indices
            @test dm.primaryData[1] == EARTH_MOON_SYSTEM.bodyData[1]
            @test dm.primaryData[2] == EARTH_MOON_SYSTEM.bodyData[2]
            # Index order is preserved in primaryData
            @test dm.primaryData[1].name == "Earth"
            @test dm.primaryData[2].name == "Moon"
            # Non-contiguous indices select correct bodies
            dm_noncont = MBD.init_dynamicsModel(EARTH_MOON_SUN_SYSTEM, [3, 1], MBD.CR3BPDynamicsModel)
            @test dm_noncont.primaryData[1].name == "Sun"
            @test dm_noncont.primaryData[2].name == "Earth"
        end

        @testset "CR3BPDynamicsModel convenience constructor" begin
            # CR3BPDynamicsModel(systemData, indices) returns CR3BPDynamicsModel
            @test MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, [1, 2]) isa MBD.CR3BPDynamicsModel
            # CR3BPDynamicsModel(systemData, indices) is consistent with init_dynamicsModel
            dm_convenience = MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, [1, 2])
            dm_initializer = MBD.init_dynamicsModel(EARTH_MOON_SYSTEM, [1, 2], MBD.CR3BPDynamicsModel)
            @test dm_convenience == dm_initializer
            # Convenience constructor propagates ArgumentError for empty indices
            @test_throws ArgumentError MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, Vector{Int64}())
        end

        @testset "CR3BP Base.show" begin
            dm = MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, [1, 2])
            # show(io, MIME, sd) does not throw
            buf = IOBuffer()
            @test_nowarn show(buf, MIME"text/plain"(), dm)
            # show output contains all body names
            out = String(take!(buf))
            @test occursin("Earth", out)
            @test occursin("Moon", out)
            # show output contains SPICE IDs
            @test occursin("399", out)
            @test occursin("301", out)
            # show output contains model type name
            @test occursin("CR3BP", out)
            # show(io, sd) delegates to MIME method without throwing
            @test_nowarn show(buf, dm)
        end

        @testset "Type hierarchy" begin
            # CR3BPDynamicsModel is a subtype of AbstractDynamicsModel
            @test MBD.CR3BPDynamicsModel <: MBD.AbstractDynamicsModel
            # Instance satisfies isa AbstractDynamicsModel
            @test MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, [1, 2]) isa MBD.AbstractDynamicsModel
        end

        @testset "AbstractDynamicsModel equality" begin
            dm1 = MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, [1, 2])
            dm2 = MBD.CR3BPDynamicsModel(EARTH_MOON_SYSTEM, [1, 2])
            # Same model compares equal
            @test dm1 == dm2
            @test isequal(dm1, dm2)
            # Different bodies are not equal
            dm_diff = MBD.CR3BPDynamicsModel(SUN_EARTH_MOON_SYSTEM, [1, 2])
            @test dm1 != dm_diff
        end
    end
end


@testset "SystemData methods" begin
    @testset "getNumPrimaries" begin
        # Returns correct count for populated SystemData
        sd1 = MBD.SystemData(["Earth", "Moon"])
        nPrimaries = getNumPrimaries(sd1)
        @test nPrimaries == 2
        @test nPrimaries isa Int64
        # Returns correct count for single-body SystemData
        sd2 = MBD.SystemData(["Earth"])
        @test getNumPrimaries(sd2) == 1
        # Returns 0 for empty SystemData
        sd_empty = MBD.SystemData(Vector{MBD.BodyData}(), Vector{String}(), Vector{Int64}())
        @test getNumPrimaries(sd_empty) == 0
        # Throws ArgumentError when names length is inconsistent
        sd_incon1 = MBD.SystemData(["Earth", "Moon"])
        push!(sd_incon1.names, "Sun")
        @test_throws ArgumentError getNumPrimaries(sd_incon1)
        # Throws ArgumentError when spiceIDs length is inconsistent
        sd_incon2 = MBD.SystemData(["Earth", "Moon"])
        push!(sd_incon2.spiceIDs, 1000)
        @test_throws ArgumentError getNumPrimaries(sd_incon2)
        # Throws ArgumentError when both names and spiceIDs lengths are inconsistent
        push!(sd_incon1.spiceIDs, 1000)
        @test_throws ArgumentError getNumPrimaries(sd_incon1)
    end

    @testset "shallowClone" begin
        # Returns a SystemData object
        sd = SystemData(["Earth", "Moon"])
        sd_copy = MBD.shallowClone(sd)
        @test sd_copy isa MBD.SystemData
        # Copy has equal field value to original
        @test isequal(sd_copy.bodyData, sd.bodyData)
        @test sd_copy.names == sd.names
        @test sd_copy.spiceIDs == sd.spiceIDs
        # Copy is a distinct SystemData object
        @test sd_copy !== sd
        @test sd_copy.bodyData !== sd.bodyData
        @test sd_copy.names !== sd.names
        @test sd_copy.spiceIDs !== sd.spiceIDs
        # Mutating copy vectors does not affect original
        push!(sd_copy.bodyData, MBD.BodyData("Sun"))
        push!(sd_copy.names, "Sun")
        push!(sd_copy.spiceIDs, 10)
        @test length(sd.bodyData) == 2
        @test length(sd.names) == 2
        @test length(sd.spiceIDs) == 2
        # Mutating original vectors does not affect copy
        push!(sd.bodyData, MBD.BodyData("Sun"), MBD.BodyData("Mars"))
        push!(sd.names, "Sun", "Mars")
        push!(sd.spiceIDs, 10, 499)
        @test length(sd_copy.bodyData) == 3
        @test length(sd_copy.names) == 3
        @test length(sd_copy.spiceIDs) == 3
        # bodyData elements are shared references (shallow, not deep)
        @test sd_copy.bodyData[1] === sd.bodyData[1]
        # Copy of empty SystemData returns empty SystemData
        sd_empty = MBD.SystemData(Vector{MBD.BodyData}(), Vector{String}(), Vector{Int64}())
        sd_emptyCopy = MBD.shallowClone(sd_empty)
        @test sd_emptyCopy isa MBD.SystemData
        @test isempty(sd_emptyCopy.bodyData)
        @test isempty(sd_emptyCopy.names)
        @test isempty(sd_emptyCopy.spiceIDs)
        # Copying is idempotent across multiple calls
        sd_copy1 = MBD.shallowClone(sd)
        sd_copy2 = MBD.shallowClone(sd)
        @test isequal(sd_copy1.bodyData, sd_copy2.bodyData)
        @test sd_copy1.names == sd_copy2.names
        @test sd_copy1.spiceIDs == sd_copy2.spiceIDs
        @test sd_copy1 !== sd_copy2
    end
end

@testset "DynamicsModel methods" begin
    @testset "Generic" begin
        struct StubDynamicsModel <: MBD.AbstractDynamicsModel
            primaryData::Vector{MBD.BodyData}
        end
        stub = StubDynamicsModel(Vector{MBD.BodyData}())

        @testset "adjustInitialConditions" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError adjustInitialConditions(stub, collect(Float64, 1:6), MBD.ARCLENGTH, MBD.FULL)
        end

        @testset "appendExtraInitialConditions" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError appendExtraInitialConditions(stub, collect(Float64, 1:6), MBD.FULL)
        end
        
        @testset "getCharLengths" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError getCharLengths(stub)
        end

        @testset "getCharMasses" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError getCharMasses(stub)
        end

        @testset "getCharTimes" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError getCharTimes(stub)
        end

        @testset "getNumPrimaries" begin
            # Returns correct count for populated dynamics model
            sd_CR3BP = MBD.SystemData(["Earth", "Moon"])
            dm_CR3BP = MBD.CR3BPDynamicsModel(sd_CR3BP, [1, 2])
            nPrimaries = getNumPrimaries(dm_CR3BP)
            @test nPrimaries == 2
            @test nPrimaries isa Int64
            # Returns 0 for empty dynamics model
            @test getNumPrimaries(stub) == 0
        end

        @testset "getMassRatios" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError getMassRatios(stub)
        end

        @testset "getStateSize" begin
            # Throws MethodError for unimplemented type
            @test_throws MethodError getStateSize(stub, MBD.FULL)
        end
    end

    @testset "CR3BP" begin
        @testset "adjustInitialConditions" begin
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            n_simple = getStateSize(dm, MBD.SIMPLE)
            n_STM = getStateSize(dm, MBD.STM)
            n_full = getStateSize(dm, MBD.FULL)
            n_arclength = getStateSize(dm, MBD.ARCLENGTH)
            q0_simple = collect(Float64, 1:n_simple)
            q0_STM = collect(Float64, 1:n_STM)
            q0_full = collect(Float64, 1:n_full)
            q0_arclength = collect(Float64, 1:n_arclength)
            # Return type is Vector{Float64}
            result = adjustInitialConditions(dm, q0_simple, MBD.SIMPLE, MBD.SIMPLE)
            @test result isa Vector{Float64}
            # Output length matches expected state size for outputEquationType
            for (q0_in, eqIn, eqOut) in [
                (q0_simple, MBD.SIMPLE, MBD.SIMPLE),
                (q0_simple, MBD.SIMPLE, MBD.STM),
                (q0_simple, MBD.SIMPLE, MBD.FULL),
                (q0_simple, MBD.SIMPLE, MBD.ARCLENGTH),
                (q0_STM,    MBD.STM,    MBD.SIMPLE),
                (q0_STM,    MBD.STM,    MBD.STM),
                (q0_STM,    MBD.STM,    MBD.FULL),
                (q0_STM,    MBD.STM,    MBD.ARCLENGTH),
                (q0_full,   MBD.FULL,   MBD.SIMPLE),
                (q0_full,   MBD.FULL,   MBD.STM),
                (q0_full,   MBD.FULL,   MBD.FULL),
            ]
                result = adjustInitialConditions(dm, q0_in, eqIn, eqOut)
                @test length(result) == getStateSize(dm, eqOut)
            end
            # Direct truncation - leading elements are preserved
            for (q0_in, eqIn, eqOut) in [
                (q0_STM,        MBD.STM,        MBD.SIMPLE),
                (q0_arclength,  MBD.ARCLENGTH,  MBD.SIMPLE),
                (q0_full,       MBD.FULL,       MBD.SIMPLE),
                (q0_full,       MBD.FULL,       MBD.STM),
            ]
                n_out = getStateSize(dm, eqOut)
                result = adjustInitialConditions(dm, q0_in, eqIn, eqOut)
                @test result[1:n_simple] == q0_in[1:n_simple]
            end
            # Indirect truncation - leading elements are preserved but next element is zero
            result_indirect = adjustInitialConditions(dm, q0_full, MBD.FULL, MBD.ARCLENGTH)
            @test result[1:n_simple] == q0_simple
            @test result[n_simple+1] == 0
            # Same-size input/output - returns leading elements unchanged
            result_same = adjustInitialConditions(dm, q0_arclength, MBD.ARCLENGTH, MBD.MOMENTUM)
            @test result[1:n_simple] == q0_arclength[1:n_simple]
            # Extension - simple states are copied to leading portion
            for (q0_in, eqIn, eqOut) in [
                (q0_simple,     MBD.SIMPLE,     MBD.STM),
                (q0_simple,     MBD.SIMPLE,     MBD.FULL),
                (q0_simple,     MBD.SIMPLE,     MBD.ARCLENGTH),
                (q0_STM,        MBD.STM,        MBD.FULL),
                (q0_arclength,  MBD.ARCLENGTH,  MBD.FULL),
            ]
                n_in = getStateSize(dm, eqIn)
                result = adjustInitialConditions(dm, q0_in, eqIn, eqOut)
                @test result[1:n_simple] == q0_in[1:n_simple]
            end
            # Extension - trailing appended elements are zero beyond STM block
            result = adjustInitialConditions(dm, q0_simple, MBD.SIMPLE, MBD.FULL)
            @test all(result[n_STM+1:end] .== 0.0)
            # Extension - STM block initialized as identity
            STM_block = reshape(result[n_simple+1:n_STM], n_simple, n_simple)
            @test STM_block ≈ LinearAlgebra.I(n_simple)
            # Extension - STM block preserved when input already contains STM
            result_STM = adjustInitialConditions(dm, q0_STM, MBD.STM, MBD.FULL)
            STM_block_STM = reshape(result_STM[n_simple+1:n_STM], n_simple, n_simple)
            expected = reshape(q0_STM[n_simple+1:n_STM], n_simple, n_simple)
            @test STM_block_STM == expected
            # Throws ArgumentError when q0 is too short for inputEquationType
            q0_short = zeros(Float64, n_simple-1)
            @test_throws ArgumentError adjustInitialConditions(dm, q0_short, MBD.SIMPLE, MBD.FULL)
            # Throws ArgumentError when q0 is too long for inputEquationType
            q0_long = zeros(Float64, n_simple+1)
            @test_throws ArgumentError adjustInitialConditions(dm, q0_long, MBD.SIMPLE, MBD.FULL)
        end

        @testset "appendExtraInitialConditions" begin
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            n_simple = getStateSize(dm, MBD.SIMPLE)
            n_STM = getStateSize(dm, MBD.STM)
            n_full = getStateSize(dm, MBD.FULL)
            n_arclength = getStateSize(dm, MBD.ARCLENGTH)
            q0_simple = collect(Float64, 1:n_simple)
            # Return type is Vector{Float64}
            @test appendExtraInitialConditions(dm, q0_simple, MBD.FULL) isa Vector{Float64}
            # Output and length matches adjustInitialConditions with SIMPLE input type
            for eqOut in [MBD.SIMPLE, MBD.STM, MBD.FULL, MBD.ARCLENGTH]
                expected = adjustInitialConditions(dm, q0_simple, MBD.SIMPLE, eqOut)
                result = appendExtraInitialConditions(dm, q0_simple, eqOut)
                @test result == expected
                @test length(result) == getStateSize(dm, eqOut)
            end
            # Propagates ArgumentError from adjustInitialConditions for wrong q0 length
            q0_bad = zeros(Float64, n_simple+1)
            @test_throws ArgumentError appendExtraInitialConditions(dm, q0_bad, MBD.FULL)
        end

        @testset "getCharLengths" begin
            # Returns secondary orbital radius for valid model
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            lstar = getCharLengths(dm)
            @test lstar == sd.bodyData[2].a
            # Return type is Float64
            @test lstar isa Float64
            # Return value is positive and finite
            @test isfinite(lstar)
            @test lstar > 0.0
            # Throws DomainError when secondary a is NaN
            bd_nan = MBD.BodyData(NaN, 1.0, 1.0, 1.0, "NaNBody", 399, 1.0, 1001, 1.0, 1.0)
            dm_nan = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_nan])
            err1 = try
                getCharLengths(dm_nan)
                nothing
            catch e
                e
            end
            @test err1 isa DomainError
            @test occursin("finite", err1.msg)
            # Throws DomainError when secondary a is non-positive
            bd_neg = MBD.BodyData(-1.0, 1.0, 1.0, 1.0, "NegBody", 399, 1.0, 1002, 1.0, 1.0)
            dm_neg = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_neg])
            err2 = try
                getCharLengths(dm_neg)
                nothing
            catch e
                e
            end
            @test err2 isa DomainError
            @test occursin("positive", err2.msg)
        end

        @testset "getCharMasses" begin
            # Returns sum of gravitational parameters divided by the gravitational constant for valid model
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            mstar = getCharMasses(dm)
            @test mstar == (sd.bodyData[1].μ+sd.bodyData[2].μ)/MBD.GRAVITY
            # Return type is Float64
            @test mstar isa Float64
            # Return value is positive and finite
            @test isfinite(mstar)
            @test mstar > 0.0
            # Throws DomainError when primary μ is NaN
            bd_nan1 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NaNPrimary", 10, 1.0, 399, NaN, 1.0)
            dm_nan1 = MBD.CR3BPDynamicsModel([bd_nan1, MBD.BodyData("Moon")])
            err1 = try
                getCharMasses(dm_nan1)
                nothing
            catch e
                e
            end
            @test err1 isa DomainError
            @test occursin("μ_1", err1.msg)
            @test occursin("finite", err1.msg)
            # Throws DomainError when secondary μ is NaN
            bd_nan2 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NaNSecondary", 399, 1.0, 1001, NaN, 1.0)
            dm_nan2 = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_nan2])
            err2 = try
                getCharMasses(dm_nan2)
                nothing
            catch e
                e
            end
            @test err2 isa DomainError
            @test occursin("μ_2", err2.msg)
            @test occursin("finite", err2.msg)
            # Throws DomainError when both μs are NaN
            dm_nan_both = MBD.CR3BPDynamicsModel([bd_nan1, bd_nan2])
            err3 = try
                getCharMasses(dm_nan_both)
                nothing
            catch e
                e
            end
            @test err3 isa DomainError
            @test occursin("μ_1", err3.msg)
            # Throws DomainError when primary μ is non-positive
            bd_neg1 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NegPrimary", 10, 1.0, 399, -1.0, 1.0)
            dm_neg1 = MBD.CR3BPDynamicsModel([bd_neg1, MBD.BodyData("Moon")])
            err4 = try
                getCharMasses(dm_neg1)
                nothing
            catch e
                e
            end
            @test err4 isa DomainError
            @test occursin("μ_1", err4.msg)
            @test occursin("positive", err4.msg)
            # Throws DomainError when secondary μ is non-positive
            bd_neg2 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NegSecondary", 399, 1.0, 1001, -1.0, 1.0)
            dm_neg2 = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_neg2])
            err5 = try
                getCharMasses(dm_neg2)
                nothing
            catch e
                e
            end
            @test err5 isa DomainError
            @test occursin("μ_2", err5.msg)
            @test occursin("positive", err5.msg)
            # Throws DomainError when both μs are non-positive
            dm_neg_both = MBD.CR3BPDynamicsModel([bd_neg1, bd_neg2])
            err6 = try
                getCharMasses(dm_neg_both)
                nothing
            catch e
                e
            end
            @test err6 isa DomainError
            @test occursin("μ_1", err6.msg)
        end

        @testset "getCharTimes" begin
            # Returns orbit period scale derived from Kepler's 3rd law for valid model
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            tstar = getCharTimes(dm)
            @test tstar == sqrt(sd.bodyData[2].a^3/(sd.bodyData[1].μ+sd.bodyData[2].μ))
            # Return type is Float64
            @test tstar isa Float64
            # Return value is positive and finite
            @test isfinite(tstar)
            @test tstar > 0.0
            # Throws DomainError when primary μ is NaN
            bd_nan1 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NaNPrimary", 10, 1.0, 399, NaN, 1.0)
            dm_nan1 = MBD.CR3BPDynamicsModel([bd_nan1, MBD.BodyData("Moon")])
            err1 = try
                getCharTimes(dm_nan1)
                nothing
            catch e
                e
            end
            @test err1 isa DomainError
            @test occursin("μ_1", err1.msg)
            @test occursin("finite", err1.msg)
            # Throws DomainError when secondary μ is NaN
            bd_nan2 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NaNSecondary", 399, 1.0, 1001, NaN, 1.0)
            dm_nan2 = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_nan2])
            err2 = try
                getCharTimes(dm_nan2)
                nothing
            catch e
                e
            end
            @test err2 isa DomainError
            @test occursin("μ_2", err2.msg)
            @test occursin("finite", err2.msg)
            # Throws DomainError when both μs are NaN
            dm_nan_both = MBD.CR3BPDynamicsModel([bd_nan1, bd_nan2])
            err3 = try
                getCharTimes(dm_nan_both)
                nothing
            catch e
                e
            end
            @test err3 isa DomainError
            @test occursin("μ_1", err3.msg)
            # Throws DomainError when primary μ is non-positive
            bd_neg1 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NegPrimary", 10, 1.0, 399, -1.0, 1.0)
            dm_neg1 = MBD.CR3BPDynamicsModel([bd_neg1, MBD.BodyData("Moon")])
            err4 = try
                getCharTimes(dm_neg1)
                nothing
            catch e
                e
            end
            @test err4 isa DomainError
            @test occursin("μ_1", err4.msg)
            @test occursin("positive", err4.msg)
            # Throws DomainError when secondary μ is non-positive
            bd_neg2 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NegSecondary", 399, 1.0, 1001, -1.0, 1.0)
            dm_neg2 = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_neg2])
            err5 = try
                getCharTimes(dm_neg2)
                nothing
            catch e
                e
            end
            @test err5 isa DomainError
            @test occursin("μ_2", err5.msg)
            @test occursin("positive", err5.msg)
            # Throws DomainError when both μs are non-positive
            dm_neg_both = MBD.CR3BPDynamicsModel([bd_neg1, bd_neg2])
            err6 = try
                getCharTimes(dm_neg_both)
                nothing
            catch e
                e
            end
            @test err6 isa DomainError
            @test occursin("μ_1", err6.msg)
        end

        @testset "getMassRatios" begin
            # Returns ratio of secondary to total mass for valid model
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            μ = getMassRatios(dm)
            @test μ == sd.bodyData[2].μ/(sd.bodyData[1].μ+sd.bodyData[2].μ)
            # Return type is Float64
            @test μ isa Float64
            # Return value is positive and finite
            @test isfinite(μ)
            @test μ > 0.0
            # Throws DomainError when primary μ is NaN
            bd_nan1 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NaNPrimary", 10, 1.0, 399, NaN, 1.0)
            dm_nan1 = MBD.CR3BPDynamicsModel([bd_nan1, MBD.BodyData("Moon")])
            err1 = try
                getMassRatios(dm_nan1)
                nothing
            catch e
                e
            end
            @test err1 isa DomainError
            @test occursin("μ_1", err1.msg)
            @test occursin("finite", err1.msg)
            # Throws DomainError when secondary μ is NaN
            bd_nan2 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NaNSecondary", 399, 1.0, 1001, NaN, 1.0)
            dm_nan2 = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_nan2])
            err2 = try
                getMassRatios(dm_nan2)
                nothing
            catch e
                e
            end
            @test err2 isa DomainError
            @test occursin("μ_2", err2.msg)
            @test occursin("finite", err2.msg)
            # Throws DomainError when both μs are NaN
            dm_nan_both = MBD.CR3BPDynamicsModel([bd_nan1, bd_nan2])
            err3 = try
                getMassRatios(dm_nan_both)
                nothing
            catch e
                e
            end
            @test err3 isa DomainError
            @test occursin("μ_1", err3.msg)
            # Throws DomainError when primary μ is non-positive
            bd_neg1 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NegPrimary", 10, 1.0, 399, -1.0, 1.0)
            dm_neg1 = MBD.CR3BPDynamicsModel([bd_neg1, MBD.BodyData("Moon")])
            err4 = try
                getMassRatios(dm_neg1)
                nothing
            catch e
                e
            end
            @test err4 isa DomainError
            @test occursin("μ_1", err4.msg)
            @test occursin("positive", err4.msg)
            # Throws DomainError when secondary μ is non-positive
            bd_neg2 = MBD.BodyData(1.0, 1.0, 1.0, 1.0, "NegSecondary", 399, 1.0, 1001, -1.0, 1.0)
            dm_neg2 = MBD.CR3BPDynamicsModel([MBD.BodyData("Earth"), bd_neg2])
            err5 = try
                getMassRatios(dm_neg2)
                nothing
            catch e
                e
            end
            @test err5 isa DomainError
            @test occursin("μ_2", err5.msg)
            @test occursin("positive", err5.msg)
            # Throws DomainError when both μs are non-positive
            dm_neg_both = MBD.CR3BPDynamicsModel([bd_neg1, bd_neg2])
            err6 = try
                getMassRatios(dm_neg_both)
                nothing
            catch e
                e
            end
            @test err6 isa DomainError
            @test occursin("μ_1", err6.msg)
        end

        @testset "getStateSize" begin
            # Returns correct state size for each supported EquationType
            sd = MBD.SystemData(["Earth", "Moon"])
            dm = MBD.CR3BPDynamicsModel(sd, [1, 2])
            expected = Dict{MBD.EquationType, Int64}(
                MBD.SIMPLE => 6,
                MBD.STM => 42,
                MBD.ARCLENGTH => 7,
                MBD.MOMENTUM => 7,
                MBD.FULL => 44
            )
            for (eqType, expectedSize) in expected
                @test getStateSize(dm, eqType) == expectedSize
            end
            # Return type is Int64 and positive for each supported EquationType
            for eqType in instances(MBD.EquationType)
                if haskey(MBD._CR3BP_state_sizes, eqType)
                    size = getStateSize(dm, eqType)
                    @test size isa Int64
                    @test size > 0
                end
            end
            # Throws ArgumentError for unsupported EquationType
            @test_throws ArgumentError getStateSize(dm, MBD.TEST)
        end
    end
end

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
