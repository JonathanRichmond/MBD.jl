"""
CR3BP manifold arc wrapper

Author: Jonathan Richmond
C: 1/12/25
"""

import MBD: CR3BPPManifoldArc

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
