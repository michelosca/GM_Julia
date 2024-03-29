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

module FunctionTerms 

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: kb 
using SharedData: r_elastic, r_diffusion
using EvaluateExpressions: ReplaceExpressionValues
using WallFlux: DensWallFluxFunction, TempWallFluxFunction
using PowerInput: PowerInputFunction

###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - GetDensRateFunction
# - GetTempRateFunction


function GetDensRateFunction(temp::Vector{Float64}, dens::Vector{Float64},
    s::Species, species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    dens_funct = 0.0

    # Check if species is meant to have a density ODE 
    if s.has_dens_eq
        s_id = s.id
        # Loop over the reaction setreaction_lists species s 
        for r in reaction_list

            # Get species index in involved_species list
            s_index = findall( x -> x == s_id, r.involved_species )
            if s_index == Int64[]
                continue
            else
                s_index = s_index[1]
            end

            # PARTICLE PRODUCTION rates
            # Terms due to particle gain/loss, e.g. recombination, ionization
            sign = r.species_balance[s_index]
            if (sign != 0) 
                if r.case == r_diffusion
                    value = sign * dens[r.reactant_species[1]] * r.K_value
                else
                    value = sign * prod(dens[r.reactant_species]) * r.K_value
                end
                dens_funct += value
                #print("   - Gain loss: ", r.id," - ", value, "\n")
            end
        end

        if s.has_wall_loss
            value = DensWallFluxFunction(s, system)
            dens_funct += value 
            #print("   - Flux loss:     ", value, "\n")
        end

        if s.has_flow_rate
            dens_funct += s.flow_rate / system.V
        end

    end
    return dens_funct
end


function GetTempRateFunction(temp::Vector{Float64}, dens::Vector{Float64},
    s::Species, species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID, t_sim::Float64)

    temp_funct = 0.0 
    if s.has_temp_eq
        # "Constants" that are used later
        Q0 = 1.5 * kb * s.dens
        Q1 = -1.5 * kb * s.temp
        s_id = s.id

        # Loop over the reaction set
        #S_abs = 0.0
        #S_elast = 0.0
        #S_inelast = 0.0
        #S_mass = 0.0
        #S_flux = 0.0
        for r in reaction_list

            # Check whether species is involved in current reaction 
            s_index = findall( x -> x == s_id, r.involved_species )
            if s_index == Int64[]
                continue
            else
                s_index = s_index[1]
            end

            # PARTICLE PRODUCTION rates
            # Terms due to particle gain/loss, e.g. recombination, ionization
            sign = r.species_balance[s_index]
            if (sign != 0)
                value = sign * prod(dens[r.reactant_species]) * r.K_value * Q1 / Q0
                temp_funct += value
                #print("   - Gain loss: ", r.id," - ", value, "\n")
                #S_mass += value
            else
                # Add energy term due to elastic collisions
                if (r.case == r_elastic)
                    # Set the neutral species mass
                    n_id = r.neutral_species_id[1] # Neutral species id
                    m_neutral = species_list[n_id].mass 
                    m_charged = s.mass
                    Q2 = -3.0 * kb * m_charged / m_neutral
                    t_neutral = temp[n_id]
                    value = Q2 * prod(dens[r.reactant_species]) *
                        r.K_value * (temp[s_id] - t_neutral) / Q0
                    temp_funct += value
                    #print("   - Elastic:   ", r.id," - ", value, "\n")
                    #S_elast += value
                end

            end

            # COLLISION INTRINSIC ENERGY GAIN/LOSS
            Er = r.E_threshold
            if !(Er == 0)
                # Check that species is part of the reacting species
                s_index = findall( x -> x == s_id, r.reactant_species)
                if s_index == Int64[]
                    continue
                else
                    s_index = s_index[1]
                end

                value = -Er * prod(dens[r.reactant_species]) * r.K_value / Q0
                temp_funct += value
                #print("   - Ethreshold: ", r.id," - ", value, "\n")
                #S_inelast += value
            end
        end

        if s.has_wall_loss
            value = TempWallFluxFunction(temp, s, species_list, system,
                system.plasma_potential, sID) / Q0
            temp_funct += value
            #print("   - Flux loss: ", value, "\n")
            #S_flux += value
        end

        if s.has_heating_mechanism
            value = PowerInputFunction(s, system, sID, t_sim) / Q0
            temp_funct += value
            #print("   - Heating :  ", value, "\n")
            #S_abs += value
        end
        #S_tot = S_abs + S_elast + S_inelast + S_mass + S_flux
        #print(S_tot,"; S_abs: ",S_abs,"; S_elast: ",S_elast,"; S_inelast: ", S_inelast,"; S_mass: ",S_mass,"; S_flux: ",S_flux,"\n")
    end
    return temp_funct
end

end