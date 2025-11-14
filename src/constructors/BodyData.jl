"""
BodyData functions

Author: Jonathan Richmond
C: 11/14/25
"""

import MBD:BodyData
import LightXML, Logging

export load_bodyData


"""
Load BodyData for a named body from an XML body data file.

Arguments
- `name::String`: The case-insensitive body name (as used by SPICE) to load.
- `fileName::String`: Path to an XML file containing `<body>` entries. When omitted
    callers often use the convenience constructor `BodyData(name::String)` which
    supplies the package's `body_data.xml`.

Returns
- `BodyData`: A fully-populated `BodyData` instance for the requested body.

Errors
- Throws `ArgumentError` if `name` or `fileName` are empty, or if the file does
    not exist.
- Throws `ErrorException` for SPICE lookup failures, XML parse errors, missing
    or empty required tags, invalid numeric fields, or when the requested body
    entry cannot be found in the XML.

Notes
- A `<parentId>` value of the literal string `"NaN"` (case-insensitive) is
    treated as "no parent" and is represented using `MBD.UNINITIALIZED_INDEX`.
- The function emits `Logging` messages at `@debug`, `@info`, `@warn`, and
    `@error` levels to aid troubleshooting.

Example
```
bd = load_bodyData("Earth", joinpath(@__DIR__, "body_data.xml"))
```
"""
function load_bodyData(name::String, fileName::String)::BodyData
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
                v = parse(Int16, parentText)
                Logging.@debug "Parsed parentId" parentSPICEID=v name=name
                v
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
