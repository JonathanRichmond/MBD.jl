"""
SystemData constructors

Author: Jonathan Richmond
C: 11/21/25
"""

import MBD:SystemData
import Logging

export init_systemData


"""
Initialize a `SystemData` container for a list of body names.

Arguments
- `names::Vector{String}`: Vector of body names to include in the system. Each
    entry must be a non-empty `String` and all entries must be unique.

Returns
- `SystemData`: A container holding the loaded `BodyData` objects, the
    original `names` vector, and the corresponding vector of `SPICEIDs`.

Errors
- Throws `ArgumentError` if `names` is empty, contains non-string or empty
    entries, or contains duplicate names.
- Throws `ErrorException` if loading any individual `BodyData(name)` fails.

Notes
- The function emits `Logging` messages at `@debug`, `@info`, `@warn`, and
    `@error` levels to aid troubleshooting.

Example
```
sys = init_systemData(["Earth", "Moon"])
```
"""
function init_systemData(names::Vector{String})::SystemData
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
