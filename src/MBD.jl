"""
Multi-body dynamics astrodynamics package

Author: Jonathan Richmond
C: 11/21/25
U: 12/12/25
"""
module MBD


import Base: show, ==
import LightXML, Logging, Printf, SPICE


const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0

include("dynamics/SystemDataMethods.jl")
include("types/BodyData.jl")
include("types/Core.jl")
include("types/SystemData.jl")
include("types/DynamicsModel.jl")
include("utilities/SPICE.jl")


export
    # Types
    BodyData,
    SystemData,
    CR3BPDynamicsModel,

    # Functions
    getIDCode


end # module MBD
