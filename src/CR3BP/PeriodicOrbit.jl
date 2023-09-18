"""
CR3BP periodic orbit wrapper

Author: Jonathan Richmond
C: 1/16/23
U: 9/10/23
"""

import LinearAlgebra
import MBD: CR3BPPeriodicOrbit

export getManifold, getStability!

"""
    getManifold(periodicOrbit, stabilitity, d)

Return stable or unstable manifold tubes

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `stability::String`: Desired manifold stability
- `d::Float64`: Stepoff distance [ndim]
"""
function getManifold(periodicOrbit::CR3BPPeriodicOrbit, stability::String, d::Float64)
    index::Int64 = (stability == "Stable") ?  argmax(abs(periodicOrbit.eigenvalues)) : argmin(abs(periodicOrbit.eigenvalues))
    eigenvector::Vector{Complex{Float64}} = periodicOrbit.eigenvectors[index]
    normEigenvector::Vector{Complex{Float64}} = eigenvector./LinearAlgebra.norm(eigenvector[1:3])
    posManifoldArc = MBD.CR3BPManifoldArc(periodicOrbit.initialCondition+d.*normEigenvector, periodicOrbit)
    negManifoldArc = MBD.CR3BPManifoldArc(periodicOrbit.initialCondition-d.*normEigenvector, periodicOrbit)

    return ([posManifoldArc], [negManifoldArc])
end

"""
    getStability!(periodicOrbit)

Return periodic orbit object with updated stability properties

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
"""
function getStability!(periodicOrbit::CR3BPPeriodicOrbit)
    E::LinearAlgebra.Eigen = LinearAlgebra.eigen(periodicOrbit.monodromy)
    periodicOrbit.eigenvalues = E.values
    periodicOrbit.eigenvectors = E.vectors
    alpha::Float64 = 2-LinearAlgebra.tr(periodicOrbit.monodromy)
    beta::Float64 = 0.5*((alpha^2)+2-LinearAlgebra.tr(periodicOrbit.monodromy^2))
    periodicOrbit.BrouckeStability = [alpha, beta]
    periodicOrbit.nu = LinearAlgebra.norm(periodicOrbit.eigenvalues, Inf)
    periodicOrbit.tau = periodicOrbit.period/log(periodicOrbit.nu)
end
