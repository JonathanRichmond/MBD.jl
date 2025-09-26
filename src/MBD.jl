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
@enum ModelType CR3BP BCR4BP EMPTY KP


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
Base.show(io::IO, bodyData::BodyData) = print(io, "BodyData for ", bodyData.name)
Base.show(io::IO, ::MIME"text/plain", bodyData::BodyData) = begin
    println(io, "BodyData")
    println(io, "\tName: ", bodyData.name, " (SPICE ID: ", bodyData.SPICEID, ")")
    println(io, "\tParent SPICE ID: ", bodyData.parentSPICEID)
    println(io, "\tGravitational parameter: ", bodyData.gravParam, " [kg^3/s^2]")
    println(io, "\tMass: ", bodyData.mass, " [kg]")
    println(io, "\tRadius: ", bodyData.bodyRadius, " [km]")
    println(io, "\tOrbit radius: ", bodyData.orbitRadius, " [km]")
    println(io, "\tInclination: ", bodyData.inc, " [rad]")
    println(io, "\tRAAN: ", bodyData.RAAN, " [rad]")
end

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
        expectedPrimaries::Dict{MBD.ModelType, Int64} = Dict{MBD.ModelType, Int64}(
            MBD.EMPTY => 0,
            MBD.KP => 1,
            MBD.CR3BP => 2,
            MBD.BCR4BP => 3
        )
        if haskey(expectedPrimaries, modelType)
            expectedNumPrimaries::Int64 = expectedPrimaries[modelType]
            if numPrimaries != expectedNumPrimaries
                Logging.@error "Model type $modelType requires $expectedNumPrimaries primary bodies, got $numPrimaries"
                throw(ArgumentError("Model type requires $expectedNumPrimaries, got $numPrimaries"))
            end
        else
            Logging.@error "Model type $modelType not found in expected primaries dictionary"
            throw(ArgumentError("Model type not found in expected primaries dictionary"))
        end

        Logging.@info "Creating $(numPrimaries+1)-body system with primaries: $(join(primaryNames, ", "))"

        primaryData::Vector{BodyData} = [BodyData(name) for name in primaryNames]
        primarySPICEIDs::Vector{Int16} = [data.SPICEID for data in primaryData]
        if numPrimaries > 1
            for j in 2:numPrimaries
                (primaryData[j].parentSPICEID != primaryData[j-1].SPICEID) && Logging.@warn "$(primaryData[j].name) does not have $(primaryData[j-1].name) as parent"
            end
        end

        return new(modelType, primaryData, collect(primaryNames), primarySPICEIDs)
    end
end
Base.:(==)(systemData1::SystemData, systemData2::SystemData) = ((systemData1.modelType == systemData2.modelType) && (systemData1.primaryData == systemData2.primaryData) && (systemData1.primaryNames == systemData2.primaryNames) && (systemData1.primarySPICEIDs == systemData2.primarySPICEIDs))
Base.show(io::IO, systemData::SystemData) = print(io, "SystemData for ", join(systemData.primaryNames, "-"), " ", length(systemData.primaryNames)+1, "-body system")
Base.show(io::IO, ::MIME"text/plain", systemData::SystemData) = begin
    println(io, "SystemData")
    println(io, "\tModel type: ", systemData.modelType)
    for (j, bodyData) in enumerate(systemData.primaryData)
        println(io, "\tPrimary ", j, ": ", bodyData.name, " (SPICE ID: ", bodyData.SPICEID, ")")
    end
end

"""
    DynamicsModel(systemData)

Dynamics model object

# Arguments
- `systemData::SystemData`: System data object
"""
mutable struct DynamicsModel
    systemData::SystemData                                              # System data object

    function DynamicsModel(systemData::SystemData)
        Logging.@info "Creating $(length(systemData.primaryNames)+1)-body dynamics model with primaries: $(join(systemData.primaryNames, ","))"

        return new(systemData)
    end
end
Base.:(==)(dynamicsModel1::DynamicsModel, dynamicsModel2::DynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)
Base.show(io::IO, dynamicsModel::DynamicsModel) = print(io, "DynamicsModel for ", join(dynamicsModel.systemData.primaryNames, "-"), " ", length(dynamicsModel.systemData.primaryNames)+1, "-body system")
Base.show(io::IO, ::MIME"text/plain", dynamicsModel::DynamicsModel) = begin
    println(io, "DynamicsModel")
    println(io, "\tSystem data: ", dynamicsModel.systemData)
end

"""
    EquationsOfMotion(dynamicsModel, equationType)

Equations of motion object

# Arguments
- `dynamicsModel::DynamicsModel`: Dynamics model object
- `equationType::EquationType`: Enumerated equation type
"""
struct EquationsOfMotion
    dynamicsModel::DynamicsModel                                        # Dynamics model object
    equationType::EquationType                                          # Enumerated equation type

    function EquationsOfMotion(dynamicsModel::DynamicsModel, equationType::EquationType)
        Logging.@info "Creating $(length(dynamicsModel.systemData.primaryNames)+1)-body $equationType EoMs with primaries: $(join(dynamicsModel.systemData.primaryNames, ","))"

        return new(dynamicsModel, equationType)
    end
end
Base.:(==)(equations1::EquationsOfMotion, equations2::EquationsOfMotion) = ((equations1.dynamicsModel == equations2.dynamicsModel) && (equations1.equationType == equations2.equationType))
Base.show(io::IO, equations::EquationsOfMotion) = print(io, "EquationsOfMotion for ", join(equations.dynamicsModel.systemData.primaryNames, "-"), " ", length(equations.dynamicsModel.systemData.primaryNames)+1, "-body system")
Base.show(io::IO, ::MIME"text/plain", equations::EquationsOfMotion) = begin
    println(io, "EquationsOfMotion")
    println(io, "\tDynamics model: ", equations.dynamicsModel)
    println(io, "\tEquation type: ", equations.equationType)
end

"""
    Integrator(integratorType)

Integrator object

# Arguments
- `integratorType::IntegratorType`: Enumerated integrator type
"""
mutable struct Integrator
    integrator                                                          # Integrator object
    integratorType::IntegratorType                                      # Enumerated integrator type

    function Integrator(integratorType::IntegratorType)
        integratorMap = Dict{IntegratorType, Any}(
            AB5     => DifferentialEquations.AB5(),
            ABM54   => DifferentialEquations.ABM54(),
            BS5     => DifferentialEquations.BS5(),
            DP5     => DifferentialEquations.DP5(),
            DP8     => DifferentialEquations.DP8(),
            VERN9   => DifferentialEquations.Vern9()
        )
        Logging.@info "Creating $integratorType integrator"

        integrator = get(integratorMap, integratorType, nothing)
        if integrator === nothing
            Logging.@error "Unsupported integrator type: $integratorType"
            throw(ArgumentError("Invalid integrator type: $integratorType"))
        end

        return new(integrator, integratorType)
    end
end
Base.:(==)(integrator1::Integrator, integrator2::Integrator) = ((integrator1.integrator == integrator2.integrator) && (integrator1.integratorType == integrator2.integratorType))
Base.show(io::IO, integrator::Integrator) = print(io, "Integrator for ", integrator.integratorType)
Base.show(io::IO, ::MIME"text/plain", integrator::Integrator) = begin
    println(io, "Integrator")
    println(io, "\tType: ", integrator.integratorType)
end

"""
    Propagator(; integrator, equationType)

Propagator object

# Arguments
- `integrator::Integrator`: Integrator object (default = Integrator(DP8))
- `equationType::EquationType`: Enumerated equation type (default = SIMPLE)
"""
mutable struct Propagator
    absTol::Float64                                                     # Absolute tolerance
    equationType::EquationType                                          # Enumerated equation type
    events::Vector{Any}                                                 # Propagation events
    integrator::Integrator                                              # Integrator object
    maxEvaluations::Int64                                               # Maximum equation evaluations
    maxStep::Int64                                                      # Maximum step size
    relTol::Float64                                                     # Relative tolerance

    function Propagator(; integrator::Integrator = Integrator(DP8), equationType::EquationType = SIMPLE)
        Logging.@info "Creating $(integrator.integratorType) propagator for $equationType EoMs"

        return new(1E-12, equationType, [], integrator, typemax(Int64), 100, 1E-12)
    end
end
Base.:(==)(propagator1::Propagator, propagator2::Propagator) = ((propagator1.equationType == propagator2.equationType) && (propagator1.integrator == propagator2.integrator))
Base.show(io::IO, propagator::Propagator) = print(io, "Propagator for ", propagator.equationType, " EoMs using ", propagator.integrator.integratorType)
Base.show(io::IO, ::MIME"text/plain", propagator::Propagator) = begin
    println(io, "Propagator")
    println(io, "\tEquation type: ", propagator.equationType)
    println(io, "\tIntegrator: ", propagator.integrator)
    println(io, "\tAbsolute tolerance: ", propagator.absTol)
    println(io, "\tRelative tolerance: ", propagator.relTol)
    println(io, "\tMaximum step: ", propagator.maxStep)
    println(io, "\tMaximum evaluations: ", propagator.maxEvaluations)
    println(io, "\tEvents: ", propagator.events)
end

"""
    Arc(dynamicsModel)

Arc object

# Arguments
- `dynamicsModel::DynamicsModel`: Dynamics model object
"""
mutable struct Arc
    dynamicsModel::DynamicsModel                                        # Dynamics model object
    states::Vector{Vector{Float64}}                                     # State vectors along arc
    times::Vector{Float64}                                              # Times along arc

    function Arc(dynamicsModel::DynamicsModel)
        Logging.@info "Creating $(dynamicsModel.systemData.modelType) arc"

        return new(dynamicsModel, [], [])
    end
end
Base.:(==)(arc1::Arc, arc2::Arc) = ((arc1.dynamicsModel == arc2.dynamicsModel) && (arc1.states == arc2.states) && (arc1.times == arc2.times))
Base.show(io::IO, arc::Arc) = print(io, "Arc for ", arc.dynamicsModel.systemData.modelType, " system")
Base.show(io::IO, ::MIME"text/plain", arc::Arc) = begin
    println(io, "Arc")
    println(io, "\tDynamics model: ", arc.dynamicsModel)
    isempty(arc.times) && return println(io, "\tNo states/times available") 
    println(io, "\tInitial state: ", arc.states[1])
    println(io, "\tInitial time: ", arc.times[1])
    println(io, "\tFinal state: ", arc.states[end])
    println(io, "\tFinal time: ", arc.times[end])
end

include("dynamics/SystemData.jl")


end # module MBD
