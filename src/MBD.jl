"""
Multi-body dynamics astrodynamics package

Author: Jonathan Richmond
C: 6/25/25
"""
module MBD

import Base: ==
import Combinatorics, DifferentialEquations, LightXML, LinearAlgebra, Logging, SPICE, StaticArrays


const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0


"""
Enumerated type for EOMs
"""
@enum EquationType ARCLENGTH FULL MOMENTUM SIMPLE STM

"""
Enumerated type for integrators
"""
@enum IntegratorType AB5 ABM54 BS5 DP5 DP8 VERN9

"""
Enumerated type for models
"""
@enum ModelType CR3BP BCR4BP TBP


"""
    BodyData(bodyName)

Body data object

# Arguments
- `bodyName::String`: Name of body
"""
mutable struct BodyData
    bodyRadius::Float64                                                 # Body mean radius [km]
    gravParam::Float64                                                  # Gravitational parameter [kg^3/s^2]
    inc::Float64                                                        # Orbit inclination relative to Ecliptic J2000 frame [rad]
    mass::Float64                                                       # Mass [kg]
    name::String                                                        # Name
    orbitRadius::Float64                                                # Mean orbit radius [km]
    parentSPICEID::Int16                                                # Parent body SPICE ID
    RAAN::Float64                                                       # Orbit Right-Ascension of Ascending Node (RAAN) relative to Ecliptic J2000 frame [rad]
    SPICEID::Int16                                                      # SPICE ID

    function BodyData(bodyName::String)
        dataFile::String = joinpath(@__DIR__, "body_data.xml")
        Logging.@info "Initializing BodyData for '$bodyName'"
        
        ID::Int16 = try
            Int16(SPICE.bods2c(bodyName))
        catch e
            Logging.@error "Could not resolve SPICE ID for body name: '$bodyName': $(sprint(showerror, e))"
            rethrow(e)
        end
        Logging.@debug "Resolved body name '$bodyName' to SPICE ID: $ID"

        try
            doc::LightXML.XMLDocument = LightXML.parse_file(dataFile)
            root::LightXML.XMLElement = LightXML.root(doc)
            nodeList::Vector{LightXML.XMLElement} = LightXML.get_elements_by_tagname(root, "body")
            
            for node in nodeList
                getValue(tag) = LightXML.content(LightXML.get_elements_by_tagname(node, tag)[1])

                SPICEID::Int16 = parse(Int16, getValue("id"))
                Logging.@debug "Examining node with SPICE ID: $SPICEID"
                
                if SPICEID == ID
                    Logging.@info "Found matching entry for SPICE ID: $ID"

                    name::String = getValue("name")
                    parentID::String = getValue("parentId")
                    gm::Float64 = parse(Float64, getValue("gm"))
                    radius::Float64 = parse(Float64, getValue("radius"))
                    circ_r::Float64 = parse(Float64, getValue("circ_r"))
                    inc::Float64 = parse(Float64, getValue("inc"))
                    RAAN::Float64 = parse(Float64, getValue("raan"))
                    parentSPICEID::Int16 = (parentID == "NaN") ? Int16(-1) : parse(Int16, parentID)
                    mass::Float64 = gm/GRAVITY

                    Logging.@debug "Parsed values: name=$name, gm=$gm [kg^3/s^2], radius=$radius [km], orbit radius=$circ_r [km], inc=$inc [rad], RAAN=$RAAN [rad]"

                    return new(radius, gm, inc, mass, name, circ_r, parentSPICEID, RAAN, SPICEID)
                end
            end

            Logging.@error "Body with SPICE ID $ID not found in XML"
            throw(ArgumentError("Body with ID $ID not found in '$dataFile'"))
        catch e
            Logging.@error "Failed to parse '$dataFile': $(sprint(showerror, e))"
            rethrow(e)
        end
    end
end
Base.:(==)(bodyData1::BodyData, bodyData2::BodyData) = (bodyData1.SPICEID == bodyData2.SPICEID)

"""
    SystemData(modelType, primaryNames)

System data object

# Arguments
- `modelType::ModelType`: Enumerated model type
- `primaryNames::Vararg{String}`: Primary names
"""
mutable struct SystemData
    modelType::ModelType                                                # Enumerated model type
    primaryData::Vector{BodyData}                                       # Primary body data objects
    primaryNames::Vector{String}                                        # Primary names
    primarySPICEIDs::Vector{Int16}                                      # Primary SPICE IDs

    function SystemData(modelType::ModelType, primaryNames::Vararg{String})
        numPrimaries::Int64 = length(primaryNames)
        Logging.@info "Creating $(numPrimaries+1)-body system with primaries: $(join(primaryNames, ", "))"

        primaryData::Vector{BodyData} = [BodyData(name) for name in primaryNames]
        primarySPICEIDs::Vector{Int16} = [data.SPICEID for data in primaryData]
        for j in 2:numPrimaries
            (primaryData[j].parentSPICEID != primaryData[j-1].SPICEID) && Logging.@warn "$(primaryData[j].name) does not have $(primaryData[j-1].name) as parent"
        end

        return new(modelType, primaryData, collect(primaryNames), primarySPICEIDs)
    end
end
Base.:(==)(systemData1::SystemData, systemData2::SystemData) = ((systemData1.modelType == systemData2.modelType) && (systemData1.primaryData == systemData2.primaryData) && (systemData1.primaryNames == systemData2.primaryNames) && (systemData1.primarySPICEIDs == systemData2.primarySPICEIDs))


end # module MBD
