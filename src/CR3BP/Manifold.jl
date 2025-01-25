"""
CR3BP manifold wrapper

Author: Jonathan Richmond
C: 1/25/25
"""

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
