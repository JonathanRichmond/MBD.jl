"""
Multi-body dynamics astrodynamics package

Author: Jonathan Richmond
C: 6/25/22
"""
module MBD

import Base: ==
import Combinatorics, DifferentialEquations, LightXML, LinearAlgebra, SPICE, StaticArrays

const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0

end # module MBD
