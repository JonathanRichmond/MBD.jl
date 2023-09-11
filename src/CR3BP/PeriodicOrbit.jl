"""
CR3BP periodic orbit wrapper

Author: Jonathan Richmond
C: 1/16/23
U: 9/10/23
"""

import LinearAlgebra
import MBD: CR3BPPeriodicOrbit

export getStability!

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
