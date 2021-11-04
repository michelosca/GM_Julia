module GenerateODEs

using SharedData: System, Species, Reaction
using SharedData: kb 
using InputBlock_Reactions: r_elastic_id, r_case_energy_sink, r_wall_loss
using WallFlux: DensWallFluxFunction, TempWallFluxFunction
using PowerInput: PowerInputFunction

###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - GenerateDensRateFunctionList
#   - GenerateDensRateFunction
#     - DensWallFluxFunction
# - GenerateTempRateFunctionList
#   - GenerateTempRateFunction
#     - TempWallFluxFunction
#     - PowerInputFunction


function GenerateDensRateFunctionList(system::System,
    species_list::Vector{Species}, reaction_list::Vector{Reaction})
    # Gathers the Density rate functions from GenerateDensRateFunction for 
    #each species in a list
    
    # Loop over all species
    dens_rate_funct_list = []
    for s in species_list
        sdens_funct_list = GenerateDensRateFunction(s, species_list,
            reaction_list, system)
        if (length(sdens_funct_list)>0)
            push!(dens_rate_funct_list, sdens_funct_list)
        end
    end
    return dens_rate_funct_list
end


function GenerateDensRateFunction(s::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction},
    system::System)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    sdens_funct_list = Function[]

    # Check if species is meant to have a density ODE 
    if (s.has_dens_eq)
        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s 
            if (r.case == r_case_energy_sink)
                continue
            end
            s_index = findall( x -> x == s.id, r.involved_species )
            if !(length(s_index) == 0)
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                # Only a density gain/loss term is added if there is a non-zero
                # particle balance
                if !(sign == 0)
                    push!(sdens_funct_list,
                        function (dens::Vector{Float64}, temp::Vector{Float64})
                            return sign * prod(dens[r.reactant_species]) *
                            r.rate_coefficient(temp)
                        end
                    )
                end
            end
        end

        if (s.has_wall_loss)
            push!(sdens_funct_list,
                function (dens::Vector{Float64}, temp::Vector{Float64})
                    return DensWallFluxFunction(dens, temp, s, species_list,
                        reaction_list, system)
                end
            )
        end
    else
        push!(sdens_funct_list, (dens, temp) -> 0)
    end
    return sdens_funct_list
end


function GenerateTempRateFunctionList(system::System,
    species_list::Vector{Species}, reaction_list::Vector{Reaction})
    # Generate temperature equations
    # Loop over all species
    temp_rate_funct_list = []
    for s in species_list
        stemp_funct_list = GenerateTempRateFunction(s, species_list,
            reaction_list, system)
        if (length(stemp_funct_list)>0)
            push!(temp_rate_funct_list, stemp_funct_list)
        end
    end
    return temp_rate_funct_list
end


function GenerateTempRateFunction(s::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

    stemp_funct_list = Function[]
    if (s.has_temp_eq)
        s_id = s.id
        # "Constants" that are used later
        Q0(dens::Vector{Float64}) = 3.0/2.0 * kb * dens[s_id]
        Q1(temp::Vector{Float64}) = -3.0/2.0 * kb * temp[s_id]

        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s_id
            s_index = findall( x -> x == s_id, r.involved_species )
            if r.id == r_wall_loss
                push!(stemp_funct_list,
                    function (dens::Vector{Float64}, temp::Vector{Float64})
                        wall_loss_rate = GetWallLossRate(dens, temp, r, system) 
                        return wall_loss_rate / Q0(dens)
                    end
                )
            elseif !(length(s_index) == 0)
                # Species s is involved in reaction r
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                if (sign != 0)
                    # Temp gain/loss term is added if there is a non-zero particle balance
                    push!(stemp_funct_list,
                        function (dens::Vector{Float64}, temp::Vector{Float64})
                            return sign * prod(dens[r.reactant_species]) *
                                r.rate_coefficient(temp) *
                                Q1(temp) / Q0(dens)
                        end
                    )
                else
                    # Add energy term due to elastic collisions
                    if (r.id == r_elastic_id)
                        # Find the neutral collision partner
                        neutral_id = r.neutral_species_id 
                        # Set the neutral species mass
                        m_neutral = species_list[neutral_id].mass 
                        m_charged = s.mass
                        Q2 = -3.0 * kb * m_charged / m_neutral
                        t_neutral = species_list[neutral_id].temp0
                        push!(stemp_funct_list,
                            function (dens::Vector{Float64}, temp::Vector{Float64})
                                e_rate = Q2 * prod(dens[r.reactant_species]) *
                                    r.rate_coefficient(temp) * (temp[s_id] - t_neutral)
                                return e_rate / Q0(dens)
                            end
                        )
                    end

                end

                # COLLISION INTRINSIC ENERGY GAIN/LOSS
                Er = r.E_threshold
                if (Er != 0) 
                    push!(stemp_funct_list,
                        function (dens::Vector{Float64}, temp::Vector{Float64})
                            e_rate = - Er * prod(dens[r.reactant_species]) *
                                r.rate_coefficient(temp)
                            return e_rate / Q0(dens)
                        end
                    )
                end
            end
        end

        if (s.has_wall_loss)
            push!(stemp_funct_list,
                function (dens::Vector{Float64}, temp::Vector{Float64})
                    flux = TempWallFluxFunction(dens, temp, s, species_list,
                        reaction_list, system)
                    return flux / Q0(dens)
                end
            )
        end

        if (s.has_heating_mechanism)
            push!(stemp_funct_list,
                function (dens::Vector{Float64}, temp::Vector{Float64})
                    ipower = PowerInputFunction(s, system) 
                    return ipower / Q0(dens)
                end
            )
        end
    else
        push!(stemp_funct_list, (dens, temp) -> 0)
    end
    return stemp_funct_list
end


end