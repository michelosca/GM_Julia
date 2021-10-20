module SharedData

struct Species
    # Species identification flags
    id::Int64          # ID for this specific species. Links to the dens, temp vectors
    species_id::Int64  # In case of ions or excited states, the ID of the
                       # neutral species. If it does not apply then = 0

    # Total rate coefficient 
    r_elastic_id::Int64

    # Physical constants
    mass::Float64
    charge::Float64

    # Dens and temp ODE flags
    has_dens_eq::Bool
    has_temp_eq::Bool
    
    # Wall flux functions
    has_wall_loss::Bool

    # Input power mechanism
    has_heating_mechanism::Bool
end

struct Reaction
    # Reaction identifier
    id::Int64
    neutral_species_id::Int64

    # Species involved in the reaction
    involved_species::Vector{Int}
    species_balance::Vector{Int}
    reactant_species::Vector{Int}
    #product_species::Vector{Int}

    # Rate coefficient function
    rate_coefficient::Function

    # Energy threshold
    E_threshold::Float64
end


# Physical constants
const kb = 1.380649e-23
const e = 1.602176634e-19     # C
const me = 9.1093837015e-31   # kg
const mp = 1.67262192369e-27  # kg
const eps0 = 8.854187817e-12

# Power input methods
const p_ccp_id = 1
const p_icp_id = 2

# Define reaction ids
const r_elastic_id = 1
const r_excitat_id = 2
const r_ionizat_id = 3
const r_recombi_id = 4

###############################################################################
# SYSTEM PARAMETERS
global A = 0.0                              # m^2
global V = 0.0                              # m^3
global l = 0.0                              # m
global drivf = 0.0                          # Hz
global drivOmega = 0.0                      # rad/s
global power_input_method = 0
global drivP = 0                            # W 
global drivV = 0                            # V
global drivI = 0                            # Amps

function SetSystemParameters(name, var, units)

    # Initialize system parameters
    lname = lowercase(name)
    if (name=="A" || lname=="area")
        global A = parse(Float64,var)
        return 0 
    elseif (name=="V" || lname=="volume")
        global V = parse(Float64, var)
        return 0
    elseif (lname=="L" || lname=="length")
        global l = parse(Float64, var) 
        return 0
    elseif (name=="f" || lname=="frequency" ||
        lname=="freq" || lname=="driving_frequency")
        if (units=="mhz")
            unit = 1.e6
        elseif (units=="ghz")
            unit = 1.e9
        elseif (units=="khz")
            unit = 1.e3
        else
            unit = 1.0
        end
        global drivf = parse(Float64,var) * unit
        global drivOmega = drivf * 2.0 * π
        return 0
    elseif (name=="P" || lname=="power_input" ||
        lname=="power")
        if (var=="from_I_and_V")
            global drivP = drivV * drivI
        else
            if (units=="kw")
                unit = 1.e3
            else
                unit = 1
            end
            global drivP = parse(Float64, var) * unit
        end
        return 0
    elseif (lname=="voltage")
        if (var=="from_I_and_P")
            global drivV = drivP / drivI
        else
            global drivV = parse(Float64, var)
        end
        return 0
    elseif (lname=="power_method" || lname=="input_power_method")
        lvar = lowcase(var)
        if (lvar == "ccp")
            p_id = p_ccp_id
        elseif (lvar == "icp")
            p_id = p_icp_id
        else
            p_id = 0
        end
        global power_input_method = p_id
        return 0
    elseif (name == "I" || lname == "current" ||
        lname == "driving_current")
        if (var=="from_P_and_V")
            global drivI = drivP / drivV
        else
            global drivI = parse(Float64,var)
        end
        return 0
    end
    # In case there is no match -> return 1
    return 1
end


# Species and reaction lists involves in the GM model
global species_list = Species[]
global reaction_list = Reaction[]

function AddReactionToList(id::Int64, invol_s::Vector{Int64},
    react_s::Vector{Int64}, balan_s::Vector{Int64}, K::Function, E::Float64,
    n_id::Int64)

    # Add react to reaction_list
    react = Reaction(id, n_id, invol_s, balan_s, react_s, K, E) 
    push!(reaction_list, react)
end


function AddSpeciesToList(id::Int64, mass::Float64, charge::Float64,
    neq_flag::Bool, Teq_flag::Bool, wl_flag::Bool, P_flag::Bool,
    n_id::Int64, r_ela_id::Int64)
    
    # Add react to reaction_list
    species = Species(id, n_id, r_ela_id, mass, charge, neq_flag,
        Teq_flag, wl_flag, P_flag) 
    push!(species_list, species)
end


# Set species ids
global s_electron_id = 0
global s_Ar_id = 0
global s_ArIon_id = 0

function SetSpeciesID(species_name, id)
    if ("e" == species_name || "electrons" == species_name)
        global s_electron_id = id
        return 0
    elseif ("Ar" == species_name)
        global s_Ar_id = id
        return 0
    elseif ("Ar+" == species_name)
        global s_ArIon_id = id
        return 0 
    end

    # In case there is no match -> return 1
    return 1
end

end