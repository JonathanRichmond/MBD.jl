"""
Multi-body dynamics astrodynamics package

Author: Jonathan LeFevre Richmond
C: 4/14/26
U: 4/30/26
"""
module MBD


import Base: show, ==
import LightXML, LinearAlgebra, Logging, Printf, SPICE


const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0


include("types/Core.jl")
include("types/BodyData.jl")
include("types/SystemData.jl")
# include("types/DynamicsModel.jl")
# include("types/EquationsOfMotion.jl")

include("dynamics/SystemDataMethods.jl")
# include("dynamics/DynamicsModelMethods.jl")
# include("dynamics/EquationsOfMotionMethods.jl")
include("utilities/spice.jl")

include("exports.jl")


end # module MBD
