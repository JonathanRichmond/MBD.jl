"""
CR3BP manifold arc wrapper

Author: Jonathan Richmond
C: 1/25/25
U: 1/27/25
"""

import MBD: CR3BPManifoldArc

export getJacobiConstant

"""
    getJacobiConstant(manifoldArc)

Return Jacobi constant

# Arguments
- `manifoldArc::CR3BPManifoldArc`: CR3BP manifold arc object
"""
function getJacobiConstant(manifoldArc::CR3BPManifoldArc)
    return getJacobiConstant(manifoldArc.periodicOrbit)
end
