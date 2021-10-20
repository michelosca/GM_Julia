module GenerateODEs

using SharedData: kb
using SharedData: Species
using SharedData: r_elastic_id
using SharedData: species_list, reaction_list
using WallFluxModels: DensWallFluxFunction, TempWallFluxFunction
using PowerInput: PowerInputFunction

function GenerateDensRateFunction(s::Species)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    sdens_funct_list = Function[]

    # Check if species is meant to have a density ODE 
    if (s.has_dens_eq)
        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s 
            s_index = findall( x -> x == s.id, r.involved_species )
            if !(length(s_index) == 0)
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                # Only a density gain/loss term is added if there is a non-zero
                # particle balance
                if !(sign == 0)
                    push!(sdens_funct_list, (dens::Vector{Float64},
                        temp::Vector{Float64}) ->
                        sign * prod(dens[r.reactant_species]) *
                        r.rate_coefficient(temp))
                end
            end
        end

        if (s.has_wall_loss)
            push!(sdens_funct_list, (dens::Vector{Float64},
                temp::Vector{Float64}) -> DensWallFluxFunction(dens, temp, s))
        end
    else
        push!(sdens_funct_list, () -> 0)
    end
    return sdens_funct_list
end


function GenerateDensRateFunctionList()
    # Gathers the Density rate functions from GenerateDensRateFunction for 
    #each species in a list
    
    # Loop over all species
    dens_rate_funct_list = []
    for s in species_list
        sdens_funct_list = GenerateDensRateFunction(s)
        if (length(sdens_funct_list)>0)
            push!(dens_rate_funct_list, sdens_funct_list)
        end
    end
    return dens_rate_funct_list
end


function GenerateTempRateFunction(s::Species)

    stemp_funct_list = Function[]
    if (s.has_temp_eq)
        s_id = s.id
        # "Constants" that are used later
        Q0(dens::Vector{Float64}) = 3.0/2.0 * kb * dens[s_id]
        Q1(temp::Vector{Float64}) = 3.0/2.0 * kb * temp[s_id]

        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s_id
            s_index = findall( x -> x == s_id, r.involved_species )
            if !(length(s_index) == 0)
                # Species s is involved in reaction r
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                if (sign != 0)
                    # Temp gain/loss term is added if there is a non-zero particle balance
                    push!(stemp_funct_list, (dens::Vector{Float64},
                        temp::Vector{Float64}) ->
                        sign * prod(dens[r.reactant_species]) *
                        r.rate_coefficient(temp) * Q1(temp) )
                else
                    # Add energy term due to elastic collisions
                    if (r.id == r_elastic_id)
                        # Find the neutral collision partner
                        neutral_id = r.neutral_species_id 
                        # Set the neutral species mass
                        m_neutral = species_list[neutral_id].mass 
                        m_charged = s.mass
                        Q2 = -3 * kb * m_charged / m_neutral
                        push!(stemp_funct_list, (dens::Vector{Float64},
                            temp::Vector{Float64}) -> Q2 *
                            prod(dens[r.reactant_species]) *
                            r.rate_coefficient(temp) *
                            (temp[s_id] - temp[neutral_id]) )
                    end

                end

                # COLLISION INTRINSIC ENERGY GAIN/LOSS
                Er = r.E_threshold
                if (Er != 0) 
                    push!(stemp_funct_list, (dens::Vector{Float64}, temp::Vector{Float64}) ->
                        - Er * prod(dens[r.reactant_species]) * r.rate_coefficient(temp) )
                end
            end
        end

        if (s.has_wall_loss)
            push!(stemp_funct_list, (dens::Vector{Float64},
                temp::Vector{Float64}) -> TempWallFluxFunction(dens, temp, s))
        end

        if (s.has_heating_mechanism)
            push!(stemp_funct_list, (dens::Vector{Float64},
                temp::Vector{Float64}) -> PowerInputFunction(dens, temp, s))
        end
    else
        push!(stemp_funct_list, () -> 0)
    end
    return stemp_funct_list
end


function GenerateTempRateFunctionList()
    # Generate temperature equations
    # Loop over all species
    temp_rate_funct_list = []
    for s in species_list
        stemp_funct_list = GenerateTempRateFunction(s)
        if (length(stemp_funct_list)>0)
            push!(temp_rate_funct_list, stemp_funct_list)
        end
    end
    return temp_rate_funct_list
end


function GatherListOfFunctions(dens::Vector{Float64}, temp::Vector{Float64},
    species::Species, funct_list::Vector{Function})
    f_out = 0

    # Collision gain/loss
    for current_f in funct_list
        print("  - Value: ", current_f(dens,temp) ,"\n")
        f_out += current_f(dens, temp)
    end

    return f_out
end

end