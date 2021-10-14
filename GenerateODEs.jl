module GenerateODEs

#include("SharedData.jl")
#using .SharedData

export GenerateDensRateFunctionList, GenerateTempRateFunctionList 


function GenerateDensRateFunction(s, reaction_list)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    sdens_funct_list = Function[]

    # Check if species is meant to have a density ODE 
    if (s.has_dens_eq)
        s_id = s.species_id
        
        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s_id
            s_index = findall( x -> x == s_id, r.involved_species )
            if !(length(s_index) == 0)
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                # Only a density gain/loss term is added if there is a non-zero particle balance
                if !(sign == 0)
                    push!(sdens_funct_list, (dens::Vector{Float64}, temp::Vector{Float64}) -> 
                        sign * prod(dens[r.reactant_species]) * r.rate_coefficient(temp))
                end
            end
        end

        # WALL LOSS rates
        # Terms due to loss with walls
        if (s.has_dens_wall_loss)
            push!(sdens_funct_list, s.dens_wall_funct)
        end
    end
    return sdens_funct_list
end


#function GenerateDensODEs(species_list::Vector{Species}, reaction_list::Vector{Reaction})
function GenerateDensRateFunctionList(species_list, reaction_list)
    # Gathers the Density rate functions from GenerateDensRateFunction for 
    #each species in a list
    
    # Loop over all species
    dens_rate_funct_list = []
    for s in species_list
        sdens_funct_list = GenerateDensRateFunction(s, reaction_list)
        if (length(sdens_funct_list)>0)
            push!(dens_rate_funct_list, sdens_funct_list)
        end
    end
    return dens_rate_funct_list
end

function GenerateTempRateFunction(s, reaction_list)

    stemp_funct_list = Function[]
    if (s.has_temp_eq)
        s_id = s.species_id
        # "Constants" that are used later
        Q0(dens::Vector{Float64}) = 3.0/2.0 * kb * dens[s.id]
        Q1(temp::Vector{Float64}) = 3.0/2.0 * kb * temp[s.id]

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
                # Only a temp gain/loss term is added if there is a non-zero particle balance
                if !(sign == 0)
                    push!(stemp_funct_list, (dens::Vector{Float64}, temp::Vector{Float64}) ->
                        sign * prod(dens[r.reactant_species]) * r.rate_coefficient(temp) * Q1(temp) )
                end

                # COLLISION INTRINSIC ENERGY GAIN/LOSS
                # In this term it is also included the energy loss/gain due to elastic scattering
                # as a function of temperature
                Er = r.E_threshold             
                push!(stemp_funct_list, (dens::Vector{Float64}, temp::Vector{Float64}) ->
                    - Er(temp) * prod(dens[r.reactant_species]) * r.rate_coefficient(temp) )
            end
        end

        # WALL LOSS rates
        # Terms due to loss with walls
        if (s.has_temp_wall_loss)
            push!(stemp_funct_list, s.temp_wall_funct)
        end

        # HEATING MECHANISM
        # Terms due to power sources 
        if (s.has_heating_mechanism)
            push!(stemp_funct_list, s.input_power_funct)
        end
    end
    return stemp_funct_list
end

function GenerateTempRateFunctionList(species_list, reaction_list)
    # Generate temperature equations
    # Loop over all species
    temp_rate_funct_list = []
    for s in species_list
        stemp_funct_list = GenerateTempRateFunction(s, reaction_list)
        if (length(stemp_funct_list)>0)
            push!(temp_rate_funct_list, stemp_funct_list)
        end
    end
    return temp_rate_funct_list
end


function gather_list_of_functions(sdens_funct_list::Vector{Function},
    dens::Vector{Float64}, temp::Vector{Float64})
    f_out = 0
    for current_f in sdens_funct_list
        f_out += current_f(dens, temp)
    end
    return f_out
end

end