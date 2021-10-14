module SharedData

export Species, Reaction
export kb, e
export c_electron_id, c_neutral_id

struct Species
    species_id::Int64
    has_dens_eq::Bool
    has_temp_eq::Bool
    
    has_dens_wall_loss::Bool
    dens_wall_funct::Function

    has_temp_wall_loss::Bool
    temp_wall_funct::Function

    has_heating_mechanism::Bool
    input_power_funct::Function
end

struct Reaction
    reaction_id::Int64
    involved_species::Vector{Int}
    reactant_species::Vector{Int}
    #product_species::Vector{Int}
    species_balance::Vector{Int}
    rate_coefficient::Function
    E_threshold::Function
end

const kb = 1.380649e-23
const e = 1.602176634e-19

# Define species ids
const c_electron_id = 1
const c_neutral_id = 2

end