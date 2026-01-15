"""
Equations of motion types

Author: Jonathan LeFevre Richmond
C: 1/15/26
"""


abstract type AbstractEquationsOfMotion end


"""
    CR3BPEquationsOfMotion <: AbstractEquationsOfMotion

Equations of motion for the Circular Restricted 3-Body Problem (CR3BP).

Fields
- `dynamicsModel::CR3BPDynamicsModel`: The underlying CR3BP dynamics model containing
  the two primary bodies (primary and secondary) and their gravitational parameters.

Construction
- Direct construction: `CR3BPEquationsOfMotion(dynamicsModel::CR3BPDynamicsModel)`
  validates the dynamics model and creates an equations of motion instance.The
  `dynamicsModel` must be a valid `CR3BPDynamicsModel` instance containing exactly two
  primary bodies (primary and secondary).

Notes
- This type encapsulates the equations of motion for the CR3BP in a rotating coordinate
  frame centered on the system's barycenter. The dynamics are governed by the
  gravitational fields of two massive bodies in circular orbits.
- The `dynamicsModel` field contains all necessary information for numerical integration
  (body masses, characteristic scales, etc.).
- Integration and numerical methods use this structure to access the underlying model.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
eom = CR3BPEquationsOfMotion(model)
println(eom)  # displays formatted equations of motion information
```
"""
struct CR3BPEquationsOfMotion <: AbstractEquationsOfMotion
    dynamicsModel::CR3BPDynamicsModel

    function CR3BPEquationsOfMotion(dynamicsModel::CR3BPDynamicsModel)
        Logging.@debug "CR3BPEquationsOfMotion constructor called" type=typeof(dynamicsModel)
        
        # Basic validation
        if dynamicsModel === nothing || !isa(dynamicsModel, CR3BPDynamicsModel)
            Logging.@error "Invalid dynamicsModel supplied to CR3BPEquationsOfMotion" type=typeof(dynamicsModel)
            throw(ArgumentError("dynamicsModel must be a CR3BPDynamicsModel instance"))
        end
        
        # CR3BP requires exactly two bodies: primary then secondary
        if length(dynamicsModel.primaryData) != 2
            Logging.@error "CR3BP requires exactly two primaries" n=length(dynamicsModel.primaryData)
            throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
        end
        
        Logging.@info "Created CR3BPEquationsOfMotion" primary=dynamicsModel.primaryData[1].name secondary=dynamicsModel.primaryData[2].name
        return new(dynamicsModel)
    end
end
Base.:(==)(EoM1::MBD.CR3BPEquationsOfMotion, EoM2::MBD.CR3BPEquationsOfMotion) = (EoM1.dynamicsModel == EoM2.dynamicsModel)
function Base.show(io::IO, ::MIME"text/plain", eom::MBD.CR3BPEquationsOfMotion)
    println(io, "CR3BPEquationsOfMotion:")
    println(io, eom.dynamicsModel)
end
function Base.show(io::IO, eom::MBD.CR3BPEquationsOfMotion)
    Base.show(io, MIME"text/plain"(), eom)
end
