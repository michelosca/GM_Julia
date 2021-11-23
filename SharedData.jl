module SharedData

# Physical constants
const kb = 1.380649e-23
const e = 1.602176634e-19     # C
const me = 9.1093837015e-31   # kg
const amu = 1.66053904020e-27  # kg
const eps0 = 8.854187817e-12
const K_to_eV = kb / e

# Power input methods
const p_ccp_id = 1
const p_icp_id = 2

# Solving sheath potential methods
const s_ohmic_power = 1
const s_flux_balance = 2
const s_flux_interpolation = 3

# Read input flags 
const c_io_error = 1

# REACTION IDs 
const r_energy_sink = -1
const r_wall_loss = -2
const r_elastic= 1
const r_excitat= 2
const r_ionizat= 3
const r_recombi= 4
const r_cx= 5

# Reaction structure
mutable struct Reaction
    # Reaction identifier
    id::Int64
    case::Int64 # Reaction species case flag
    neutral_species_id::Vector{Int64}

    # Species involved in the reaction
    involved_species::Vector{Int}
    species_balance::Vector{Int}
    reactant_species::Vector{Int}

    # Rate coefficient function
    rate_coefficient

    # Energy threshold
    E_threshold::Float64


    Reaction() = new()
end


# Species structure
mutable struct Species
    # Species identification flags
    id::Int64          # ID for this specific species. Links to the dens, temp vectors
    species_id::Int64  # In case of ions or excited states, the ID of the
                       # neutral species. If it does not apply then = 0

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

    # Density and temperature of the species
    dens::Float64
    temp::Float64
    pressure::Float64

    # Other features
    reaction_list::Vector{Reaction}
    mfp::Float64
    v_thermal::Float64 # Thermal speed
    v_Bohm::Float64    # Bohm speed
    D::Float64         # Gudmundsson parameter
    Lambda::Float64
    h_R::Float64
    h_L::Float64
    gamma::Float64     # Sticking coefficient
    n_sheath::Float64
    flux::Float64

    # Output features
    name::String
    Species() = new()
end

# System structure
mutable struct System
    A::Float64                              # system area, m^2
    V::Float64                              # system volume, m^3
    l::Float64                              # system length, m
    radius::Float64                         # system radius, m

    power_input_method::Int64               # driving power method, flag
    Vsheath_solving_method::Int64
    drivf::Float64                          # driving frequency, Hz
    drivOmega::Float64                      # driving frequency, rad/s
    drivP::Float64                          # driving power, W 
    drivV::Float64                          # driving voltage, V
    drivI::Float64                          # criving current, Amps

    t_end::Float64                          # simulation time, seconds

    electrode_area::Float64                 # m^2

    System() = new()
end

struct Output
    output_flag_list::Vector{Int64}
end


mutable struct SpeciesID
    current_id::Int64

    electron::Int64

    Ar::Int64
    Ar_Ion::Int64
    Ar_Exc::Int64

    O::Int64
    O_negIon::Int64
    O_Ion::Int64
    O_3p::Int64
    O_1d::Int64


    O2::Int64
    O2_Ion::Int64
    O2_a1Ag::Int64

    SpeciesID() = new()
end

end