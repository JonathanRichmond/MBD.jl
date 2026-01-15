"""
System data types

Author: Jonathan LeFevre Richmond
C: 12/12/25
"""


"""
    init_systemData(names::Vector{String}) -> SystemData

Initialize a multi-body system by loading body data for each named body.

Arguments
- `names::Vector{String}`: Non-empty vector of canonical body names (e.g.,
  ["Earth", "Moon"]). Each name must be a non-empty string and will be used
  to load the corresponding `BodyData` via `BodyData(name)`. Duplicate names
  are not allowed.

Returns
- `SystemData`: A populated `SystemData` struct containing vectors of loaded
  `BodyData`, the original names, and corresponding SPICE IDs.

Errors
- Throws `ArgumentError` if:
  - `names` is empty
  - Any element in `names` is not a string or is empty/whitespace-only
  - `names` contains duplicate entries
- Throws `ErrorException` if loading `BodyData` for any name fails (e.g., name
  cannot be resolved via SPICE, or body data is missing from the XML file).

Notes
- The function emits `Logging` messages at `@debug`, `@info`, and `@error`
  levels to aid troubleshooting.
- Bodies are loaded in the order specified in `names`, and the resulting
  `bodyData`, `names`, and `SPICEIDs` vectors maintain this ordering.

Example
```
sys = init_systemData(["Earth", "Moon"])
```
"""
function init_systemData(names::Vector{String})::MBD.SystemData
    Logging.@debug "init_systemData called" n_names=length(names)

    # Validate input collection
    if isempty(names)
        Logging.@error "Empty names vector supplied to init_systemData"
        throw(ArgumentError("`names` must be a non-empty Vector{String}"))
    end
    for (idx, name) in enumerate(names)
        if !isa(name, String) || isempty(strip(name))
            Logging.@error "Invalid name entry in names vector" index=idx value=name
            throw(ArgumentError("names[$idx] must be a non-empty String"))
        end
    end

    # Check for duplicates which would likely indicate a user error
    if length(unique(names)) != length(names)
        Logging.@error "Duplicate entries found in names" names=names
        throw(ArgumentError("`names` contains duplicate entries"))
    end

    # Build vectors by pushing so we don't keep uninitialized entries
    bodyData = Vector{MBD.BodyData}()
    SPICEIDs = Vector{Int16}()
    for name in names
        Logging.@debug "Loading BodyData" name=name
        data = try
            MBD.BodyData(name)
        catch err
            Logging.@error "Failed to load BodyData" name=name err=err
            throw(ErrorException("failed to load BodyData for '$(name)': $(err)"))
        end
        push!(bodyData, data)
        push!(SPICEIDs, data.SPICEID)
        Logging.@info "Loaded BodyData" name=name SPICEID=data.SPICEID
    end

    Logging.@info "Initialized SystemData" n_bodies=length(bodyData)
    return SystemData(bodyData, names, SPICEIDs)
end


"""
    SystemData

Container type representing a multi-body system composed of loaded celestial bodies.

Fields
- `bodyData::Vector{BodyData}`: Vector of `BodyData` objects, one for each body
  in the system. Order corresponds to `names` and `SPICEIDs`.
- `names::Vector{String}`: The original list of canonical body names used to
  construct the system (e.g., ["Earth", "Moon"]).
- `SPICEIDs::Vector{Int16}`: The SPICE integer identifiers for each body, in the
  same order as `names` and `bodyData`.

Construction
- Direct construction: `SystemData(bodyData, names, SPICEIDs)`
- Convenience constructor: `SystemData(names::Vector{String})` calls
  `init_systemData(names)` to validate inputs and load each `BodyData`.

Notes
- All three vectors maintain consistent ordering: element `i` in each vector
  corresponds to the same body.
- Use the convenience constructor for typical usage to ensure proper validation
  and initialization.

Example
```
sys = SystemData(["Earth", "Moon"])
println(sys)  # displays formatted system information
```
"""
struct SystemData
    bodyData::Vector{MBD.BodyData}
    names::Vector{String}
    SPICEIDs::Vector{Int16}
end
SystemData(names::Vector{String}) = init_systemData(names)
Base.:(==)(sys1::MBD.SystemData, sys2::MBD.SystemData) = (sys1.names == sys2.names) && (sys1.bodyData == sys2.bodyData) && (sys1.SPICEIDs == sys2.SPICEIDs)
function Base.show(io::IO, ::MIME"text/plain", sys::MBD.SystemData)
    println(io, "SystemData: ", length(sys.bodyData), " bodies")
    println(io, "  Names and SPICE IDs:")
    for (idx, bd) in enumerate(sys.bodyData)
        Printf.@printf(io, "    %2d) %-19s  SPICEID=%6d  parent=%6d\n", idx, sys.names[idx], sys.SPICEIDs[idx], bd.parentSPICEID)
    end
end
function Base.show(io::IO, sys::MBD.SystemData)
    Base.show(io, MIME"text/plain"(), sys)
end
