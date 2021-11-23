module FunctionTerms 

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: kb 
using InputBlock_Reactions: r_elastic, r_wall_loss, r_energy_sink
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
    s::Species, reaction_list::Vector{Reaction}, system::System, sID::SpeciesID)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    dens_funct = 0.0

    # Check if species is meant to have a density ODE 
    if (s.has_dens_eq)
        # Loop over the reaction setreaction_lists species s 
        for r in reaction_list

            # Check if species s is involved in current reaction
            s_index = findall( x -> x == s.id, r.involved_species )
            if s_index == Int64[] 
                continue
            elseif r.case == r_energy_sink
                continue
            else
                s_index = s_index[1]
            end

            # PARTICLE PRODUCTION rates
            # Terms due to particle gain/loss, e.g. recombination, ionization
            sign = r.species_balance[s_index]
            if !(sign == 0) 
                if r.case == r_wall_loss
                    K = r.rate_coefficient(temp, s, system, sID) 
                else
                    K = r.rate_coefficient(temp, sID) 
                end
                value = sign * prod(dens[r.reactant_species]) * K
                dens_funct += value
                #print("   - Gain loss: ", r.id," - ", value, "\n")
            end
        end

        if (s.has_wall_loss)
            value = DensWallFluxFunction(s, system)
            dens_funct += value 
            #print("   - Flux loss:     ", value, "\n")
        end
    end
    return dens_funct
end


function GetTempRateFunction(temp::Vector{Float64}, dens::Vector{Float64},
    s::Species, species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, V_sheath::Float64, sID::SpeciesID)

    temp_funct = 0.0 
    if (s.has_temp_eq)
        # "Constants" that are used later
        Q0 = 3.0/2.0 * kb * s.dens
        Q1 = -3.0/2.0 * kb * s.temp
        s_id = s.id

        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s_id
            s_index = findall( x -> x == s_id, r.involved_species )
            if s_index == Int64[] 
                # Species is not involved in reaction -> move to next reaction
                continue
            else
                s_index = s_index[1]
            end

            # PARTICLE PRODUCTION rates
            # Terms due to particle gain/loss, e.g. recombination, ionization
            sign = r.species_balance[s_index]
            if (sign != 0)
                if r.case == r_wall_loss
                    K = r.rate_coefficient(temp, s, system, sID) 
                else
                    K = r.rate_coefficient(temp, sID) 
                end
                value = sign * prod(dens[r.reactant_species]) * K * Q1 / Q0
                temp_funct += value
                #print("   - Gain loss: ", r.id," - ", value, "\n")
            else
                # Add energy term due to elastic collisions
                if (r.case == r_elastic)
                    # Set the neutral species mass
                    n_id = r.neutral_species_id[1] # Neutral species id
                    m_neutral = species_list[n_id].mass 
                    m_charged = s.mass
                    Q2 = -3.0 * kb * m_charged / m_neutral
                    t_neutral = temp[n_id]
                    K = r.rate_coefficient(temp, sID) 
                    value = Q2 * prod(dens[r.reactant_species]) *
                        K * (temp[s_id] - t_neutral) / Q0
                    temp_funct += value
                    #print("   - Elastic:   ", r.id," - ", value, "\n")
                end

            end

            # COLLISION INTRINSIC ENERGY GAIN/LOSS
            # Check that species is part of the reacting species
            s_index = findall( x -> x == s_id, r.involved_species )
            if s_index == Int64[] 
                continue
            else
                s_index = s_index[1]
            end
            Er = r.E_threshold
            if (Er != 0) 
                K = r.rate_coefficient(temp, sID) 
                value = -Er * prod(dens[r.reactant_species]) * K / Q0
                temp_funct += value
                #print("   - Ethreshold: ", r.id," - ", value, "\n")
            end
        end

        if (s.has_wall_loss)
            value = TempWallFluxFunction(temp, s, species_list, system, V_sheath, sID) / Q0
            temp_funct += value
            #print("   - Flux loss: ", value, "\n")
        end

        if (s.has_heating_mechanism)
            value = PowerInputFunction(s, system, sID) / Q0
            temp_funct += value
            #print("   - Heating :  ", value, "\n")
        end
    end
    return temp_funct
end

end