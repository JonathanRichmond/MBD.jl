"""
SPICE utility functions

Author: Jonathan Richmond
C: 11/14/25
"""

import Logging, SPICE

export getIDCode


"""
Resolve a SPICE integer ID for a given body name.

Arguments
- `name::String`: The body name (as accepted by SPICE) to resolve.

Returns
- The integer SPICE ID corresponding to `name` as returned by
  `SPICE.bods2c`.

Errors
- Throws `ArgumentError` if `name` is empty.
- Throws `ErrorException` when the underlying `SPICE.bods2c` call fails or
  returns no code.

Logging
- Emits `@error` logs for invalid input or SPICE failures to aid debugging.

Example
```
id = getIDCode("Earth")
```
"""
function getIDCode_func(name::String)
    if isempty(strip(name))
        Logging.@error "Empty name supplied to getIDCode"
        throw(ArgumentError("`name` must be a non-empty string"))
    end
    code = try
        SPICE.bods2c(name)
    catch err
        Logging.@error "SPICE.bods2c failed" name=name err=err
        throw(ErrorException("SPICE lookup failed for name '$(name)': $(err)"))
    end
    if code === nothing
        Logging.@error "SPICE.bods2c returned no code" name=name
        throw(ErrorException("SPICE lookup returned no code for name '$(name)'") )
    end

    return code
end
const _getIDCode_func = Ref{Function}(getIDCode_func)
function getIDCode(name::String)
    code = _getIDCode_func[](name)
    if code === nothing
        Logging.@error "SPICE.bods2c returned no code" name=name
        throw(ErrorException("SPICE lookup returned no code for name '$(name)'") )
    end
    return code
end
