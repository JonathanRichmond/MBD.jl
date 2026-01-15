"""
SystemData methods

Author: Jonathan LeFevre Richmond
C: 12/22/25
"""


"""
    getNumPrimaries(sys::SystemData) -> Int64

Return the number of primary bodies contained in a `SystemData` instance.

Arguments
- `sys::SystemData`: Valid `SystemData` with internally consistent vectors.

Returns
- `Int64`: Count of entries in `sys.bodyData`.

Errors
- Throws `ArgumentError` if `sys` is not a `SystemData` instance.
- Throws `ErrorException` if `bodyData`, `names`, and `SPICEIDs` lengths differ.

Notes
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
sys = SystemData(["Earth", "Moon"])
n = getNumPrimaries(sys)
```
"""
function getNumPrimaries(sys::SystemData)
    Logging.@debug "getNumPrimaries called"

    if sys === nothing || !isa(sys, SystemData)
        Logging.@error "Invalid SystemData supplied to getNumPrimaries" sys=sys
        throw(ArgumentError("sys must be a SystemData instance"))
    end

    n::Int64 = length(sys.bodyData)
    if (length(sys.names) != n) || (length(sys.SPICEIDs) != n)
        Logging.@error "SystemData vectors length mismatch" bodyData_len=n names_len=length(sys.names) SPICEIDs_len=length(sys.SPICEIDs)
        throw(ErrorException("SystemData fields are inconsistent in length"))
    end

    Logging.@info "Computed number of primaries" count=n
    return n
end

"""
    shallowClone(sys::SystemData) -> SystemData

Create a shallow clone of a `SystemData` instance.

Arguments
- `sys::SystemData`: Valid `SystemData` to clone.

Returns
- `SystemData`: A new `SystemData` with copied data vectors and original names.

Errors
- Throws `ArgumentError` if `sys` is not a `SystemData` instance.
- Throws `ErrorException` if `bodyData`, `names`, and `SPICEIDs` lengths differ.

Behavior
- Copies `bodyData` and `SPICEIDs` vectors; reuses `names` vector reference.

Notes
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
sys = SystemData(["Earth", "Moon"])
clone = shallowClone(sys)
```
"""
function shallowClone(sys::SystemData)
    Logging.@debug "shallowClone called"

    if sys === nothing || !isa(sys, SystemData)
        Logging.@error "Invalid SystemData supplied to shallowClone" sys=sys
        throw(ArgumentError("sys must be a SystemData instance"))
    end

    n::Int64 = length(sys.bodyData)
    if (length(sys.names) != n) || (length(sys.SPICEIDs) != n)
        Logging.@error "SystemData vectors length mismatch" bodyData_len=n names_len=length(sys.names) SPICEIDs_len=length(sys.SPICEIDs)
        throw(ErrorException("SystemData fields are inconsistent in length"))
    end
    clone = SystemData(copy(sys.bodyData), sys.names, copy(sys.SPICEIDs))

    Logging.@info "Created shallow clone of SystemData" n_bodies=length(clone.bodyData)
    return clone
end
