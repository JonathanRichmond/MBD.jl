"""
Multi-body dynamics astrodynamics package

Author: Jonathan Richmond
C: 11/14/25
"""
module MBD


import Base: ==
import Printf


const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0


"""
Container type describing the physical and orbital properties of a body.

Fields
- `a::Float64` : Mean (circular) orbital radius (typically semimajor axis or
    equivalent circular radius) in the units used by the data file.
- `e::Float64` : Orbital eccentricity.
- `i::Float64` : Orbital inclination (radians).
- `m::Float64` : Mass of the body expressed consistently with `MBD.GRAVITY`.
- `name::String` : Canonical body name (as used with SPICE and the data file).
- `parentSPICEID::Int16` : SPICE ID of the parent body, or `MBD.UNINITIALIZED_INDEX`
    when the body has no parent.
- `r::Float64` : Physical radius of the body.
- `SPICEID::Int16` : SPICE integer identifier for this body.
- `μ::Float64` : Standard gravitational parameter (GM) for the body.
- `Ω::Float64` : Right ascension of the ascending node (RAAN), in radians.

Notes
- Instances are normally constructed via `BodyData(name::String)`, which calls
    `load_bodyData(name, joinpath(@__DIR__, "body_data.xml"))` to populate fields
    from the packaged XML data file.

Example
```
bd = BodyData("Earth")
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
BodyData(name::String) = load_bodyData(name, joinpath(@__DIR__, "body_data.xml"))
Base.:(==)(data1::BodyData, data2::BodyData) = (data1.SPICEID == data2.SPICEID) && (data1.a == data2.a) && (data1.e == data2.e) && (data1.i == data2.i) && (data1.m == data2.m) && (data1.name == data2.name) && (data1.parentSPICEID == data2.parentSPICEID) && (data1.r == data2.r) && (data1.μ == data2.μ) && (data1.Ω == data2.Ω)
function Base.show(io::IO, ::MIME"text/plain", data::BodyData)
    println(io, "BodyData: ", data.name)
    println(io, "  SPICEID: ", data.SPICEID, "   Parent SPICEID: ", data.parentSPICEID)
    Printf.@printf(io, "  a: %0.6g km   e: %0.6g   i: %0.6g rad\n", data.a, data.e, data.i)
    Printf.@printf(io, "  radius: %0.6g km   μ (GM): %0.6g km^3/s^2   mass: %0.6g kg\n", data.r, data.μ, data.m)
    Printf.@printf(io, "  Ω (RAAN): %0.6g rad\n", data.Ω)
end
function Base.show(io::IO, data::BodyData)
    Base.show(io, MIME"text/plain"(), data)
end


include("constructors/BodyData.jl")
include("utilities/SPICE.jl")


end # module MBD
