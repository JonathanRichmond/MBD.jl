"""
Multi-Body Dynamics astrodynamics package

Author: Jonathan Richmond
C: 9/1/22
U: 7/29/23
"""
module MBD

import Base: ==
import LightXML

const GRAVITY = 6.67384e-20

"""
Enumerated type for EOMs
"""
@enum EquationType FULL SIMPLE STM

"""
Abstract type for dynamics models
"""
abstract type AbstractDynamicsModel end

"""
Abstract type for EOMs
"""
abstract type AbstractEquationsOfMotion end

"""
Abstract type for system objects
"""
abstract type AbstractSystemData end

"""
    BodyName(name)

Body name object

# Arguments
- `name::String`: Name of body
"""
struct BodyName
    name::String                    # Name

    function BodyName(name)
        return new(name)
    end
end

"""
    BodyData(name)

Body object

# Arguments
- `name::String`: Name of body
"""
mutable struct BodyData
    bodyRadius::Float64             # Body mean radius [km]
    gravParam::Float64              # Gravitational parameter [kg^3/s^2]
    inc::Float64                    # Orbit inclination relative to Ecliptic J2000 frame [rad]
    mass::Float64                   # Mass [kg]
    name::String                    # Name
    orbitRadius::Float64            # Orbit mean radius [km]
    parentSpiceID::Int64            # Parent body SPICE ID
    RAAN::Float64                   # Orbit Right-Ascension of Ascending Node relative to Ecliptic J2000 frame [rad]
    spiceID::Int64                  # SPICE ID

    function BodyData(name::String)
        this = new()

        dataFile::String = joinpath(@__DIR__,"body_data.xml")
        bodyName = BodyName(name)
        id::Int64 = getIDCode(bodyName)
        try
            doc::LightXML.XMLDocument = LightXML.parse_file(dataFile)
            root::LightXML.XMLElement = LightXML.root(doc)
            nodeList::Vector{LightXML.XMLElement} = LightXML.get_elements_by_tagname(root, "body")
            for node::LightXML.XMLElement in nodeList
                if LightXML.is_elementnode(node)
                    elID::String = LightXML.content(LightXML.get_elements_by_tagname(node, "id")[1])
                    this.spiceID = parse(Int64, elID)
                    if this.spiceID == id
                        this.name = LightXML.content(LightXML.get_elements_by_tagname(node, "name")[1])
                        parentIDStr::String = LightXML.content(LightXML.get_elements_by_tagname(node, "parentId")[1])
                        (parentIDStr == "NaN") || (this.parentSpiceID = parse(Int64, parentIDStr))
                        this.gravParam = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "gm")[1]))
                        this.mass = this.gravParam/GRAVITY
                        this.bodyRadius = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "radius")[1]))
                        this.orbitRadius = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "circ_r")[1]))
                        this.inc = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "inc")[1]))
                        this.RAAN = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "raan")[1]))
                        break
                    end
                end
            end
            isdefined(this, :name) || throw(ArgumentError("Body $id not found in 'body_data.xml'"))
        catch
            throw(ErrorException("Could not parse $dataFile"))
        end

        return this
    end
end

"""
    CR3BPSystemData(p1, p2)

CR3BP system object

# Arguments
- `p1::String`: Name of first primary
- `p2::String`: Name of second primary
"""
mutable struct CR3BPSystemData <: AbstractSystemData
    charLength::Float64             # Characteristic length [km]
    charMass::Float64               # Characteristic mass [kg]
    charTime::Float64               # Characteritic time [s]
    numPrimaries::Int64             # Number of primaries that must exist in this system
    params::Vector{Float64}         # Other system parameters
    primaryNames::Vector{String}    # Primary names
    primarySpiceIDs::Vector{Int64}  # Primary SPICE IDs

    function CR3BPSystemData(p1::String, p2::String)
        this = new()

        p1Data = BodyData(p1)
        p2Data = BodyData(p2)
        this.numPrimaries = 2
        this.primaryNames = [p1Data.name, p2Data.name]
        this.primarySpiceIDs = [p1Data.spiceID, p2Data.spiceID]
        (p2Data.parentSpiceID == this.primarySpiceIDs[1]) || throw(ArgumentError("First primary must be parent of second primary"))
        totalGravParam::Float64 = p1Data.gravParam+p2Data.gravParam
        this.charLength = p2Data.orbitRadius
        this.charMass = totalGravParam/GRAVITY
        this.charTime = sqrt(this.charLength^3/totalGravParam)
        this.params = [p2Data.gravParam/totalGravParam]
        
        return this
    end
end
Base.:(==)(systemData1::CR3BPSystemData, systemData2::CR3BPSystemData) = (systemData1.primarySpiceIDs == systemData2.primarySpiceIDs)

"""
    CR3BPDynamicsModel(systemData)

CR3BP dynamics model object

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system object
"""
struct CR3BPDynamicsModel <: AbstractDynamicsModel
    systemData::CR3BPSystemData     # CR3BP system object

    function CR3BPDynamicsModel(systemData::CR3BPSystemData)
        return new(systemData)
    end
end
Base.:(==)(dynamicsModel1::CR3BPDynamicsModel, dynamicsModel2::CR3BPDynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)

"""
    CR3BPEquationsOfMotion(equationType, dynamicsModel)

CR3BP EOM object

# Arguments
- `equationType::EquationType`: EOM type
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct CR3BPEquationsOfMotion <: AbstractEquationsOfMotion
    dim::Int64                      # State vector dimension
    equationType::EquationType      # EOM type
    mu::Float64                     # CR3BP system mass ratio

    function CR3BPEquationsOfMotion(equationType::EquationType, dynamicsModel::CR3BPDynamicsModel)
        return new(getStateSize(dynamicsModel, equationType), equationType, getMassRatio(dynamicsModel.systemData))
    end
end

include("CR3BP/DynamicsModel.jl")
include("CR3BP/EquationsOfMotion.jl")
include("CR3BP/SystemData.jl")
include("spice/BodyName.jl")

end # module MBD
