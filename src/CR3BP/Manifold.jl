"""
CR3BP manifold wrapper

Author: Jonathan Richmond
C: 1/23/25
"""

import MBD: CR3BPPManifold

export getJacobiConstant

"""
    getJacobiConstant(manifold)

Return Jacobi constant

# Arguments
- `manifold::CR3BPManifold`: CR3BP manifold object
"""
function getJacobiConstant(manifold::CR3BPManifold)
    return getJacobiConstant(manifold.periodicOrbit)
end
