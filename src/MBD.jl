"""
Multi-body dynamics astrodynamics package

Author: Jonathan Richmond
C: 11/21/25
U: 12/22/25
"""
module MBD


import Base: show, ==
import LightXML, Logging, Printf, SPICE


const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0


include("types/Core.jl")
include("types/BodyData.jl")
include("types/SystemData.jl")
include("types/DynamicsModel.jl")

include("dynamics/SystemDataMethods.jl")
include("dynamics/DynamicsModelMethods.jl")
include("utilities/SPICE.jl")

include("exports.jl")


end # module MBD
