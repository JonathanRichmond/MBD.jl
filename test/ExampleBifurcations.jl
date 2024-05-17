"""
Bifurcation functions

Author: Jonathan Richmond
C: 5/16/24
"""

using MBD

export detectBifurcations!

"""
    detectBifurcations!(targeter, orbitFamily)

Return orbit family object with updated bifurcations

# Arguments
- `targeter::AbstractTargeter`: CR3BP orbit targeter object
- `orbitFamily::CR3BPOrbitFamily`: CR3BP orbit family object
"""
function detectBifurcations!(targeter::MBD.AbstractTargeter, orbitFamily::MBD.CR3BPOrbitFamily)
    alpha1::Float64 = orbitFamily.familyMembers[1].BrouckeStability[1]
    beta1::Float64 = orbitFamily.familyMembers[1].BrouckeStability[2]
    tangentSign::Int64 = sign(2*alpha1+beta1+2)
    tangentCount::Int64 = 0
    for b::Int64 in 2:length(orbitFamily.familyMembers)
        alpha::Float64 = orbitFamily.familyMembers[b].BrouckeStability[1]
        beta::Float64 = orbitFamily.familyMembers[b].BrouckeStability[2]
        if sign(2*alpha+beta+2) != tangentSign
            tangentSign *= -1
            tangentCount += 1
            println("\nFound tangent bifurcation!:")
            bifurcationOrbit::MBD.CR3BPPeriodicOrbit = getTangentBifurcationOrbit(targeter, orbitFamily, b-1, b)
            insert!(orbitFamily.familyMembers, b, bifurcationOrbit)
            eigenSort!(orbitFamily)
            push!(orbitFamily.bifurcations, MBD.Bifurcation(orbitFamily, bifurcationOrbit, b, MBD.TANGENT, tangentCount))
        end
    end
end
