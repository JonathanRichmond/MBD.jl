"""
Body data types

Author: Jonathan Richmond
C: 12/12/25
"""


"""
    load_bodyData(name::String, fileName::String) -> BodyData

Load physical and orbital properties for a celestial body from an XML data file.

Arguments
- `name::String`: The canonical name of the body (e.g., "Earth", "Moon"). This name
  is resolved to a SPICE ID using `getIDCode(name)` and must match a `<body>` entry
  in the XML file with the corresponding `<id>` element.
- `fileName::String`: Absolute path to the XML data file containing body definitions.
  Each `<body>` element must include tags: `<id>`, `<circ_r>`, `<ecc>`, `<inc>`,
  `<parentId>`, `<radius>`, `<gm>`, and `<raan>`.

Returns
- `BodyData`: A populated `BodyData` struct with the body's orbital elements, mass,
  radius, SPICE IDs, and gravitational parameter.

Errors
- Throws `ArgumentError` if:
  - `name` or `fileName` is empty or contains only whitespace
  - The specified file does not exist
- Throws `ErrorException` if:
  - SPICE ID resolution fails for the given `name`
  - XML parsing fails
  - Required XML elements are missing or contain invalid/empty values
  - No `<body>` element with matching SPICE ID is found in the file

Notes
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.
- The `<parentId>` field may be "NaN" (case-insensitive) to indicate no parent;
  in this case `parentSPICEID` is set to `MBD.UNINITIALIZED_INDEX`.
- Mass is computed as `μ / MBD.GRAVITY` where `μ` is read from `<gm>`.

Example
```
bodyData = load_bodyData("Earth", "/path/to/body_data.xml")
```
"""
function load_bodyData(name::String, fileName::String)::MBD.BodyData
    Logging.@debug "load_bodyData called" name=name fileName=fileName

    # Validate inputs
    if isempty(strip(name))
        Logging.@error "Empty name provided to load_bodyData" fileName=fileName
        throw(ArgumentError("`name` must be a non-empty string"))
    end
    if isempty(strip(fileName))
        Logging.@error "Empty fileName provided to load_bodyData" name=name
        throw(ArgumentError("`fileName` must be a non-empty string"))
    end

    # Resolve SPICE ID and validate
    SPICEID::Int16 = try
        id = Int16(getIDCode(name))
        Logging.@debug "Resolved SPICE ID" name=name SPICEID=id
        id
    catch err
        Logging.@error "Failed to resolve SPICE ID" name=name err=err
        throw(ErrorException("Failed to resolve SPICE ID for name '$(name)': $(err)"))
    end

    # Check file exists
     if !isfile(fileName)
        Logging.@error "Body data file not found" fileName=fileName
        throw(ArgumentError("Body data file not found: $(fileName)"))
    end

    # Parse XML document
    doc::LightXML.XMLDocument = try
        d = LightXML.parse_file(fileName)
        Logging.@debug "Parsed XML file" fileName=fileName
        d
    catch err
        Logging.@error "Failed to parse XML file" fileName=fileName err=err
        throw(ErrorException("Failed to parse XML file '$(fileName)': $(err)"))
    end
    root = LightXML.root(doc)
    root === nothing && begin
        Logging.@error "XML document has no root element" fileName=fileName
        throw(ErrorException("XML document '$(fileName)' has no root element"))
    end
    bodyList = LightXML.get_elements_by_tagname(root, "body")
    (bodyList === nothing) || (length(bodyList) == 0) && begin
        Logging.@warn "No <body> elements found in XML" fileName=fileName
        throw(ErrorException("No <body> elements found in '$(fileName)'"))
    end

    # Iterate bodies and find matching SPICE ID
    for body in bodyList
        if !LightXML.is_elementnode(body)
            Logging.@debug "Skipping non-element node in bodyList"
            continue
        end

        # Helper to fetch single element text and error if missing
        getTagText = function(element::LightXML.XMLElement, tag::String)
            nodes = LightXML.get_elements_by_tagname(element, tag)
            if (nodes === nothing) || (length(nodes) == 0)
                Logging.@error "Missing tag in body element" tag=tag fileName=fileName
                throw(ErrorException("Missing <$(tag)> for body in '$(fileName)'") )
            end
            text = LightXML.content(nodes[1])
            if (text === nothing) || isempty(strip(text))
                Logging.@error "Empty tag text in body element" tag=tag fileName=fileName
                throw(ErrorException("Empty <$(tag)> for body in '$(fileName)'") )
            end

            return text
        end

        # Read and validate ID
        idText = try
            getTagText(body, "id")
        catch err
            rethrow()
        end
        id::Int64 = try
            parse(Int64, strip(idText))
        catch err
            Logging.@error "Failed to parse <id> value" idText=idText fileName=fileName err=err
            throw(ErrorException("Invalid integer value for <id>: '$(idText)'") )
        end
        if id != SPICEID
            Logging.@debug "SPICE ID does not match; skipping body" found_id=id expected=SPICEID
            continue
        end

        # Found the matching body; read required fields with checks
        a::Float64 = try
            val = parse(Float64, strip(getTagText(body, "circ_r")))
            Logging.@debug "Parsed circ_r" value=val name=name
            val
        catch err
            Logging.@error "Invalid <circ_r>" name=name err=err
            throw(ErrorException("Invalid or missing <circ_r> for body '$(name)': $(err)"))
        end
        e::Float64 = try
            val = parse(Float64, strip(getTagText(body, "ecc")))
            Logging.@debug "Parsed ecc" value=val name=name
            val
        catch err
            Logging.@error "Invalid <ecc>" name=name err=err
            throw(ErrorException("Invalid or missing <ecc> for body '$(name)': $(err)"))
        end
        i::Float64 = try
            val = parse(Float64, strip(getTagText(body, "inc")))
            Logging.@debug "Parsed inc" value=val name=name
            val
        catch err
            Logging.@error "Invalid <inc>" name=name err=err
            throw(ErrorException("Invalid or missing <inc> for body '$(name)': $(err)"))
        end
        parentText = strip(getTagText(body, "parentId"))
        if lowercase(parentText) == "nan"
            Logging.@debug "parentId is NaN; using UNINITIALIZED_INDEX" name=name parentText=parentText
            parentSPICEID::Int16 = Int16(MBD.UNINITIALIZED_INDEX)
        else
            parentSPICEID = try
                val = parse(Int16, parentText)
                Logging.@debug "Parsed parentId" parentSPICEID=v name=name
                val
            catch err
                Logging.@error "Invalid <parentId> value" parentText=parentText name=name err=err
                throw(ErrorException("Invalid <parentId> for body '$(name)': '$(parentText)'"))
            end
        end
        r::Float64 = try
            val = parse(Float64, strip(getTagText(body, "radius")))
            Logging.@debug "Parsed radius" value=val name=name
            val
        catch err
            Logging.@error "Invalid <radius>" name=name err=err
            throw(ErrorException("Invalid or missing <radius> for body '$(name)': $(err)"))
        end
        μ::Float64 = try
            val = parse(Float64, strip(getTagText(body, "gm")))
            Logging.@debug "Parsed gm" value=val name=name
            val
        catch err
            Logging.@error "Invalid <gm>" name=name err=err
            throw(ErrorException("Invalid or missing <gm> for body '$(name)': $(err)"))
        end
        Ω::Float64 = try
            val = parse(Float64, strip(getTagText(body, "raan")))
            Logging.@debug "Parsed raan" value=val name=name
            val
        catch err
            Logging.@error "Invalid <raan>" name=name err=err
            throw(ErrorException("Invalid or missing <raan> for body '$(name)': $(err)"))
        end
        m::Float64 = μ/MBD.GRAVITY

        Logging.@info "Loaded body data" name=name SPICEID=SPICEID parentSPICEID=parentSPICEID
        return BodyData(a, e, i, m, name, parentSPICEID, r, SPICEID, μ, Ω)
    end

    # If we fall through, no matching body was found
    throw(ErrorException("No body with SPICE ID $(SPICEID) (name='$(name)') found in '$(fileName)'"))
end


"""
    BodyData

Container type describing the physical and orbital properties of a celestial body.

Fields
- `a::Float64`: Mean (circular) orbital radius, typically the semimajor axis or
  equivalent circular radius, in kilometers.
- `e::Float64`: Orbital eccentricity (dimensionless).
- `i::Float64`: Orbital inclination in radians.
- `m::Float64`: Mass of the body in kilograms, computed as `μ / MBD.GRAVITY`.
- `name::String`: Canonical body name as used with SPICE and the XML data file
  (e.g., "Earth", "Moon").
- `parentSPICEID::Int16`: SPICE integer ID of the parent body. Set to
  `MBD.UNINITIALIZED_INDEX` if the body has no parent (e.g., the Sun).
- `r::Float64`: Physical radius of the body in kilometers.
- `SPICEID::Int16`: SPICE integer identifier uniquely identifying this body.
- `μ::Float64`: Standard gravitational parameter (GM) for the body, in km³/s².
- `Ω::Float64`: Right ascension of the ascending node (RAAN) in radians.

Construction
- Direct construction: `BodyData(a, e, i, m, name, parentSPICEID, r, SPICEID, μ, Ω)`
- Convenience constructor: `BodyData(name::String)` loads data from the packaged
  XML file `body_data.xml` via `load_bodyData`.

Example
```
bd = BodyData("Earth")
println(bd)  # displays formatted body information
```
"""
struct BodyData
    a::Float64
    e::Float64
    i::Float64
    m::Float64
    name::String
    parentSPICEID::Int16
    r::Float64
    SPICEID::Int16
    μ::Float64
    Ω::Float64
end
BodyData(name::String) = load_bodyData(name, joinpath(@__DIR__, "../body_data.xml"))
Base.:(==)(data1::MBD.BodyData, data2::MBD.BodyData) = (data1.SPICEID == data2.SPICEID) && (data1.a == data2.a) && (data1.e == data2.e) && (data1.i == data2.i) && (data1.m == data2.m) && (data1.name == data2.name) && (data1.parentSPICEID == data2.parentSPICEID) && (data1.r == data2.r) && (data1.μ == data2.μ) && (data1.Ω == data2.Ω)
function Base.show(io::IO, ::MIME"text/plain", data::MBD.BodyData)
    println(io, "BodyData: ", data.name)
    println(io, "  SPICEID: ", data.SPICEID, "   Parent SPICEID: ", data.parentSPICEID)
    Printf.@printf(io, "  a: %0.6g km   e: %0.6g   i: %0.6g rad\n", data.a, data.e, data.i)
    Printf.@printf(io, "  radius: %0.6g km   μ (GM): %0.6g km^3/s^2   mass: %0.6g kg\n", data.r, data.μ, data.m)
    Printf.@printf(io, "  Ω (RAAN): %0.6g rad\n", data.Ω)
end
function Base.show(io::IO, data::MBD.BodyData)
    Base.show(io, MIME"text/plain"(), data)
end
