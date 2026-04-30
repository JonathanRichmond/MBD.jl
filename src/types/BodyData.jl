"""
Body data types

Author: Jonathan LeFevre Richmond
C: 4/14/26
U: 4/30/26
"""


"""
    load_bodyData(name::String, fileName::Sting) -> BodyData

Load the body data for a given body name from a given XML file

Arguments
- `name::String`: Body name (as accepted by SPICE)
- `fileName::String`: XML file path containing body definitions

Returns
- `BodyData`: BodyData object

Errors
- Throws `ArgumentError` if `name` is empty
- Throws `SystemError` if `fileName` does not exist
- Throws `KeyError` if no XML elements exist, an XML element tag does not
    exist, or the SPICE ID does not exist in `fileName`

Logging
- Emits `@error` logs for thrown errors
- Emits `@warn` logs if an XML element tag does not exist or when an
    unparseable SPICE ID is skipped
- Emits `@info` logs when body data is loaded
- Emits `@debug` logs when function is entered, if `name` is cached, when
    `fileName` is parsed, when a SPICE ID is retrieved, when a body list is
    scanned, when a body is matched, or when `fileName` is freed

Example
```
bodyData = load_bodyData("Earth", "path/to/body_data.xml")
```
"""
function load_bodyData(name::String, fileName::String)::BodyData
    Logging.@debug "Entered load_bodyData" name fileName
    
    # Input validation
    if isempty(strip(name))
        Logging.@error "Body name must be a non-empty string" name
        throw(ArgumentError("Body name cannot be empty or whitespace"))
    end
    if !isfile(fileName)
        Logging.@error "XML file not found"
        throw(SystemError("File not found: \"$fileName\""))
    end

    # Normalize name once to avoid redundant allocations
    normalized = String(strip(name))

    # Cache lookup (avoids redundant file parsing)
    if haskey(_body_cache, normalized)
        Logging.@debug "Cache hit for body data" name=normalized fileName
        return _body_cache[normalized]
    end

    Logging.@debug "Cache miss - parsing XML" name=normalized fileName

    # SPICE ID lookup
    local spiceID
    try
        spiceID = getIDCode(normalized)
    catch e
        Logging.@error "Failed to retrieve SPICE ID" name=normalized exception=(e, catch_backtrace())
        rethrow()
    end
    Logging.@debug "Retrieved SPICE ID" name=normalized spiceID

    # Parse XML file
    doc = nothing
    try
        doc = LightXML.parse_file(fileName)
    catch e
        Logging.@error "Failed to parse XML file" fileName exception=(e, catch_backtrace())
        rethrow()
    end

    # Helper function to extract text content from a child tag
    function getTagText(element::LightXML.XMLElement, tag::String)::String
        nodes = LightXML.get_elements_by_tagname(element, tag)
        if isempty(nodes)
            Logging.@warn "XML tag not found"
            throw(KeyError("Tag <$tag> not found in XML element"))
        end

        return strip(LightXML.content(nodes[1]))
    end
    
    # Parse XML data
    local result::BodyData
    try
        root = LightXML.root(doc)
        bodyList = LightXML.get_elements_by_tagname(root, "body")
        
        if isempty(bodyList)
            Logging.@error "No <body> elements found in XML file" fileName
            throw(KeyError("No <body> elements found in \"$fileName\"."))
        end

        Logging.@debug "Scanning body list" total=length(bodyList) target=normalized

        found = false
        for body in bodyList
            # Parse ID first and skip non-matching bodies immediately
            idText = try
                getTagText(body, "id")
            catch e
                Logging.@warn "Skipping <body> with missing <id>" exception=(e, catch_backtrace())
                continue
            end
            id = tryparse(Int64, idText)
            if isnothing(id)
                Logging.@warn "Skipping <body> with unparseable <id>" idText
                continue
            end
            id != spiceID && continue
            Logging.@debug "Matched body by SPICE ID" id spiceID

            # Parse remaining fields only for the matched body
            local a::Float64, e::Float64, i::Float64, parentSpiceID::Int64, r::Float64, μ::Float64, Ω::Float64
            try
                a = parse(Float64, getTagText(body, "circ_r"))
                e = parse(Float64, getTagText(body, "ecc"))
                i = parse(Float64, getTagText(body, "inc"))
                r = parse(Float64, getTagText(body, "radius"))
                μ = parse(Float64, getTagText(body, "gm"))
                Ω = parse(Float64, getTagText(body, "raan"))
            catch err
                Logging.@error "Failed to parse numeric field for body" name=normalized id exception=(err, catch_backtrace())
                rethrow()
            end
            parentText = getTagText(body, "parentId")
            if lowercase(parentText) == "nan"
                parentSpiceID = UNINITIALIZED_INDEX
            else
                try
                    parentSpiceID = parse(Int64, parentText)
                catch e
                    Logging.@error "Failed to parse <parentId> for body" name=normalized parentText exception=(e, catch_backtrace())
                    rethrow()
                end
            end
            m = μ/GRAVITY

            result = BodyData(a, e, i, m, normalized, parentSpiceID, r, spiceID, μ, Ω)
            found = true
            break
        end

        if !found
            Logging.@error "No matching body found in XML for SPICE ID" name=normalized spiceID fileName
            throw(KeyError("No body with SPICE ID $spiceID (\"$normalized\") found in \"$fileName\"."))
        end
    finally
        # Free the XML document to prevent memory leaks
        LightXML.free(doc)
        Logging.@debug "Freed XML document" fileName
    end

    # Store result in cache
    _body_cache[normalized] = result
    Logging.@info "Loaded body data" name=normalized spiceID fileName

    return result
end


"""
    BodyData(name::String)

Fields
- `a::Float64`: Mean (circular) orbital radius, typically semimajor axis [km]
- `e::Float64`: Mean orbital eccentricity [ndim]
- `i::Float64`: Mean orbital inclination [rad]
- `m::Float64`: Mass [kg]
- `name::String`: Body name (as accepted by SPICE)
- `parentSPICEID::Int64`: Parent body SPICE ID (set to 0 if no parent)
- `r::Float64`: Mean body radius [km]
- `SPICEID::Int64`: SPICE ID
- `μ::Float64`: Standard gravitational parameter [km^3/s^2]
- `Ω::Float64`: Mean right ascension of ascending node (RAAN) [rad]

Construction
- Direct: `BodyData(a, e, i, m, name, parentSPICEID, r, SPICEID, μ, Ω)`
- Convenience: `BodyData(name)` loads data from the packaged XML file
    `body_data.xml` via `load_bodyData`

Example
```
bodyData = BodyData("Earth")
println(bodyData)
```
"""
struct BodyData
    a::Float64
    e::Float64
    i::Float64
    m::Float64
    name::String
    parentSpiceID::Int64
    r::Float64
    spiceID::Int64
    μ::Float64
    Ω::Float64
end
BodyData(name::String) = load_bodyData(name, joinpath(@__DIR__, "../body_data.xml"))


# Cache for memoization
const _body_cache = Dict{String, BodyData}()


# Base functions
Base.:(==)(bodyData1::BodyData, bodyData2::BodyData) = (bodyData1.spiceID == bodyData2.spiceID) && (bodyData1.a == bodyData2.a) && (bodyData1.e == bodyData2.e) && (bodyData1.i == bodyData2.i) && (bodyData1.m == bodyData2.m) && (bodyData1.name == bodyData2.name) && (bodyData1.parentSpiceID == bodyData2.parentSpiceID) && (bodyData1.r == bodyData2.r) && (bodyData1.μ == bodyData2.μ) && (bodyData1.Ω == bodyData2.Ω)
Base.isequal(bodyData1::BodyData, bodyData2::BodyData) = isequal(bodyData1.spiceID, bodyData2.spiceID) && isequal(bodyData1.a, bodyData2.a) && isequal(bodyData1.e, bodyData2.e) && isequal(bodyData1.i, bodyData2.i) && isequal(bodyData1.m, bodyData2.m) && isequal(bodyData1.name, bodyData2.name) && isequal(bodyData1.parentSpiceID, bodyData2.parentSpiceID) && isequal(bodyData1.r, bodyData2.r) && isequal(bodyData1.μ, bodyData2.μ) && isequal(bodyData1.Ω, bodyData2.Ω)
function Base.show(io::IO, ::MIME"text/plain", bodyData::BodyData)
    println(io, "BodyData: ", bodyData.name)
    println(io, "  SPICE ID: ", bodyData.spiceID, "   Parent SPICE ID: ", bodyData.parentSpiceID)
    Printf.@printf(io, "  a: %0.6g km   e: %0.6g   i: %0.6g rad\n", bodyData.a, bodyData.e, bodyData.i)
    Printf.@printf(io, "  radius: %0.6g km   μ (GM): %0.6g km^3/s^2   mass: %0.6g kg\n", bodyData.r, bodyData.μ, bodyData.m)
    Printf.@printf(io, "  Ω (RAAN): %0.6g rad\n", bodyData.Ω)
end
function Base.show(io::IO, bodyData::BodyData)
    Base.show(io, MIME"text/plain"(), bodyData)
end
