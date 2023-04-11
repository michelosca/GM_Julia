# Copyright (C) 2021 Michel Osca Engelbrecht
#
# This file is part of GM Julia.
#
# GM Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GM Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module SharedData

using DataFrames

# Physical constants
const kb = 1.380649e-23
const e = 1.602176634e-19     # C
const me = 9.1093837015e-31   # kg
const amu = 1.66053904020e-27  # kg
const eps0 = 8.854187817e-12
const K_to_eV = kb / e

# Read input flags 
const c_io_error = 1

# REACTION IDs 
const r_elastic = 100
const r_diffusion = 101
const r_lower_threshold = 103
const r_emission_rate = 104
const r_extended = 105

# OUTPUT constants
const o_scale_lin = 200
const o_scale_log = 201
const o_single_run = 202
const o_pL = 203
const o_dens= 204
const o_temp = 205
const o_power = 206
const o_pressure = 207
const o_pressure_percent = 208
const o_frequency = 209
const o_duty_ratio = 210
const o_total_pressure = 211

# Species ID referring to ions/atoms
const heavy_species_id = 300

# INPUT DATA BLOCKS
const b_system = 301
const b_species = 302
const b_reactions = 303
const b_output = 304
const b_constants = 305

# h factor ID 
const h_classical = 401     # https://onlinelibrary.wiley.com/doi/pdf/10.1002/ppap.201600138
const h_Gudmundsson = 402  # https://doi.org/10.1088/0022-3727/33/11/311 , https://doi.org/10.1088/0963-0252/16/2/025
const h_Monahan = 403      # https://doi.org/10.1088/0963-0252/17/4/045003

# Solving sheath potential methods
const s_ohmic_power = 501
const s_flux_balance = 502
const s_flux_interpolation = 503

# Time dependent power waveforms
const p_constant = 600
const p_square = 601 


# Reaction structure
mutable struct Reaction
    # Reaction identifier
    name::String
    id::Int64
    case::Int64 # Reaction species case flag
    neutral_species_id::Vector{Int64}

    # Species involved in the reaction
    involved_species::Vector{Int}
    species_balance::Vector{Int}
    reactant_species::Vector{Int}
    product_species::Vector{Int}

    # Rate coefficient function
    rate_coefficient::Union{Float64, Expr, Function}
    K_value::Float64

    # Energy threshold
    E_threshold::Float64

    # Radiation emission
    self_absorption::Bool
    g_high::Float64 # statistical weight of the higher state (reactant)
    g_low::Float64  # statistical weight of the lower state (product)
    g_high_total::Float64 # statistical weight of the higher state (reactant)
    g_low_total::Float64  # statistical weight of the lower state (product)
    wavelength::Float64


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

    # Flags
    has_dens_eq::Bool
    has_temp_eq::Bool
    has_wall_loss::Bool
    has_heating_mechanism::Bool
    has_flow_rate::Bool

    # Density and temperature of the species
    dens::Float64
    temp::Float64
    pressure::Float64

    # Other features
    reaction_list::Vector{Reaction}
    mfp::Float64       # mean free path
    v_thermal::Float64 # Thermal speed
    v_Bohm::Float64    # Bohm speed

    # Flux losses
    n_sheath::Float64  # number density at sheath entrance
    flux::Float64      # particle flux [particles/s]

    # Diffusion losses
    D::Float64         # Diffusion coefficient
    h_R::Float64
    h_L::Float64
    gamma::Float64     # Sticking coefficient

    # External in/out flow
    flow_rate::Float64

    # Node parameters
    in_nodes::Vector{Tuple{Int64, String, Float64}}
    out_nodes::Vector{Tuple{Int64, String, Float64}}

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

    h_id::Int64                             # h factor ID (determines sheath density) 
    Vsheath_solving_method::Int64
    drivf::Float64                          # driving frequency, Hz
    drivOmega::Float64                      # driving frequency, rad/s
    drivP::Float64                          # driving power, W 
    P_absorbed::Float64                     # power absorbed per unit volume, W / m^3
    P_shape::Int64
    P_duty_ratio::Float64
    P_start::Float64
    dt_start::Float64
    on_slope::Float64
    off_slope::Float64

    plasma_potential::Float64
    total_pressure::Float64
    T_e_min::Float64
    T_e_max::Float64

    t_end::Float64                          # simulation time, seconds
    errcode::Int64

    alpha::Float64
    Lambda::Float64

    prerun::Bool
    folder::String

    log_file::String
    System() = new()
end


mutable struct OutputBlock

    case::Vector{Int64}
    species_id::Vector{Int64}
    scale::Vector{Int64}
    n_parameters::Int64
    label::String
    
    name::Vector{String}
    x::Vector{Float64}
    x_min::Vector{Float64}
    x_max::Vector{Float64}
    x_steps::Vector{Int64}

    n_data_frame::DataFrame
    T_data_frame::DataFrame
    K_data_frame::DataFrame
    param_data_frame::DataFrame

    first_dump::Bool
    OutputBlock() = new()
end


mutable struct SpeciesID

    current_id::Int64

    electron::Int64

    Ar::Int64
    Ar_Ion::Int64
    Ar_m::Int64
    Ar_r::Int64
    Ar_4p::Int64

    O::Int64
    O_negIon::Int64
    O_Ion::Int64
    O_1d::Int64
    O_1s::Int64
    O_3s::Int64
    O_5s::Int64
    O_3p::Int64
    O_5p::Int64


    O2::Int64
    O2_v::Int64
    O2_Ion::Int64
    O2_negIon::Int64
    O2_a1Ag::Int64
    O2_b1Su::Int64
    O2_a1Ag_v::Int64
    O2_b1Su_v::Int64

    O3::Int64
    O3_v::Int64
    O3_Ion::Int64
    O3_negIon::Int64
    O4::Int64
    O4_Ion::Int64
    O4_negIon::Int64

    SpeciesID() = new()
end 

end