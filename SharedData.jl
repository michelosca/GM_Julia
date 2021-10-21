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
    rate_coefficient::Expr

    # Energy threshold
    E_threshold::Float64
end


# Physical constants
const kb = 1.380649e-23
const e = 1.602176634e-19     # C
const me = 9.1093837015e-31   # kg
const mp = 1.67262192369e-27  # kg
const eps0 = 8.854187817e-12
const K_to_eV = kb / e

# Power input methods
const p_ccp_id = 1
const p_icp_id = 2

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

function SetSystemParameters(name, var, unit)

    errcode = 1
    # Initialize system parameters
    lname = lowercase(name)
    if (name=="A" || lname=="area")
        global A = parse(Float64,var)
        errcode = 0
    elseif (lname=="v" || lname=="volume" || lname=="vol")
        global V = parse(Float64, var)
        errcode = 0
    elseif (lname=="l" || lname=="length")
        global l = parse(Float64, var) 
        errcode = 0
    elseif (name=="f" || lname=="frequency" ||
        lname=="freq" || lname=="driving_frequency")
        global drivf = parse(Float64,var) * unit
        global drivOmega = drivf * 2.0 * π
        errcode = 0
    elseif (name=="P" || lname=="power_input" ||
        lname=="power")
        if (var=="from_I_and_V")
            global drivP = drivV * drivI
        else
            global drivP = parse(Float64, var) * unit
        end
        errcode = 0
    elseif (lname=="voltage")
        if (var=="from_I_and_P")
            global drivV = drivP / drivI
        else
            global drivV = parse(Float64, var)
        end
        errcode = 0
    elseif (lname=="power_method" || lname=="input_power_method" ||
            lname=="power_input_method")
        lvar = lowercase(var)
        if (lvar == "ccp")
            p_id = p_ccp_id
            errcode = 0
        elseif (lvar == "icp")
            p_id = p_icp_id
            errcode = 05.0e-39 * Te_eV ^ 4.5
        else
            p_id = 0
            errcode = 1 
        end
        global power_input_method = p_id
    elseif (name == "I" || lname == "current" ||
        lname == "driving_current")
        if (var=="from_P_and_V")
            global drivI = drivP / drivV
        else
            global drivI = parse(Float64,var)
        end
        errcode = 0
    end
    return errcode 
end


# Species and reaction lists involves in the GM model
global species_list = Species[]
global reaction_list = Reaction[]

function AddReactionToList(id::Int64, invol_s::Vector{Int64},
    react_s::Vector{Int64}, balan_s::Vector{Int64}, K::Expr, E::Float64,
    n_id::Int64)

    errcode = 0
    try
        # Add react to reaction_list
        react = Reaction(id, n_id, invol_s, balan_s, react_s, K, E) 
        push!(reaction_list, react)
    catch
        print("***ERROR*** While attaching reaction\n")
        errcode = 1
    end
    return errcode
end


function AddSpeciesToList(id::Int64, mass::Float64, charge::Float64,
    neq_flag::Bool, Teq_flag::Bool, wl_flag::Bool, P_flag::Bool,
    n_id::Int64, r_ela_id::Int64)
    
    errcode = 0
    try
        # Add react to reaction_list
        species = Species(id, n_id, r_ela_id, mass, charge, neq_flag,
            Teq_flag, wl_flag, P_flag) 
        push!(species_list, species)
    catch
        print("***ERROR*** While attaching species\n")
        errcode = 1
    end

    return errcode
end


# Set species ids
global s_electron_id = 0
global s_Ar_id = 0
global s_ArIon_id = 0
global s_ArExc_id = 0

function SetSpeciesID(species_name, id)

    errcode = 0

    if ("e" == species_name || "electrons" == species_name)
        global s_electron_id = id
    elseif ("Ar" == species_name)
        global s_Ar_id = id
    elseif ("Ar+" == species_name)
        global s_ArIon_id = id
    elseif ("Ar*" == species_name || "Ar excited" == species_name)
        global s_ArExc_id = id
    else
        errcode = 1
    end
    return errcode 
end

# Define reaction ids
const r_elastic_id = 1
const r_excitat_id = 2
const r_ionizat_id = 3
const r_recombi_id = 4
const r_cx_id = 5

end