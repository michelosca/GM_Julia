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

# Read input flags 
const c_io_error = 1


# Species structure
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

    # Initial conditions
    dens0::Float64
    temp0::Float64
end

# Reaction structure
struct Reaction
    # Reaction identifier
    id::Int64
    neutral_species_id::Int64

    # Species involved in the reaction
    involved_species::Vector{Int}
    species_balance::Vector{Int}
    reactant_species::Vector{Int}

    # Rate coefficient function
    rate_coefficient

    # Energy threshold
    E_threshold::Float64
end


# System structure
struct System
    A::Float64                              # system area, m^2
    V::Float64                              # system volume, m^3
    l::Float64                              # system length, m

    power_input_method::Int64               # driving power method, flag
    drivf::Float64                          # driving frequency, Hz
    drivOmega::Float64                      # driving frequency, rad/s
    drivP::Float64                          # driving power, W 
    drivV::Float64                          # driving voltage, V
    drivI::Float64                          # criving current, Amps
end

end