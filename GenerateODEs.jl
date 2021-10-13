module GenerateODEs

include("SharedData.jl")
using .SharedData

export GenerateDensODEs, GenerateTempODEs

#function GenerateDensODEs(species_list::Vector{Species}, reaction_list::Vector{Reaction})
function GenerateDensODEs(species_list, reaction_list)
    # Generate density equations
    # Loop over all species
    dens_eq_list = []
    for s in species_list
        print("Species ", s.species_id,"\n")

        # PARTICLE PRODUCTION
        # First, add (if exist) the gain/loss term due to particle production
        if (s.has_dens_eq)
            funct_list = Function[]
            s_id = s.species_id
            # Loop over the reaction set
            r_counter = 0
            for r in reaction_list
                print("  - Reaction ",r.reaction_id,"\n")
                # Check whether reaction r involves species s_id
                s_index = findall( x -> x == s_id, r.reacting_species )
                print("    - Is current species involved ",s_index,"\n")
                if !(length(s_index) == 0)
                    r_id = r.reaction_id
                    s_index = s_index[1]
                    print("    - Species ", s_id, " is part of reaction ", r_id,"\n")
                    sign = r.species_balance[s_index]
                    print("    - The species balance is: ",sign,"\n")
                    # Only a density gain/loss term is added if there is a non-zero particle balance
                    if !(sign == 0)
                        # Generate gain/loss density rate function
                        f(args...) = sign * args[3] * args[4] * r.rate_coefficient(args...)
                        print("    - Function added\n")
                        push!(funct_list, f)
                        r_counter += 1
                    end
                end
            end
            # All the terms for species s are collected.
            # Now, add them to the equation list as a single function 
            f_gainloss(args...) = sum(funct_list[i](args...) for i=1:r_counter)
        else
            # If species has no density equation then the function returns zero
            f_gainloss(args...) = 0
        end

        # WALL LOSS rates
        # Second, find the WALL LOSS for this specific species
        if (s.has_dens_wall_loss)
            f_wallloss(args...) = s.dens_wall_loss(args...)
        else
            f_wallloss(args...) = 0
        end

        # This is f in the density equation, i.e. in d(n)/dt = f(Te, Tg, ne, nAr,...)
        f_dens(args...) = f_gainloss(args...) + f_wallloss(args...)

        # Add f_dens to the list of equations
        push!(dens_eq_list, f_dens)

    end
    return dens_eq_list
end


function GenerateTempODEs(species_list, reaction_list)
    # Generate temperature equations
    # Loop over all species
    temp_eq_list = []
    for s in species_list

        # First, add (if exist) the gain/loss term due to particle production
        Q0(args...) = 3.0/2.0 * kb * args[3]
        Q1(args...) = -3.0/2.0 * kb * args[1]
        if (s.has_temp_eq)
            funct_list = Function[]
            s_id = s.species_id
            # Loop over the reaction set
            r_counter = 0
            for r in reaction_list
                # Check whether reaction r involves species s_id
                s_index = findall( x -> x == s_id, r.reacting_species )
                if !(length(s_index) == 0)
                    # Species s is involved in reaction r
                    s_index = s_index[1]

                    # PARTICLE GAIN/LOSS
                    # Energy gain/loss dur to particle gain/loss only if non-zero particle balance
                    sign = r.species_balance[s_index]
                    if !(sign == 0)
                        # Generate gain/loss temperature rate function
                        f_part(args...) = sign * args[3] * args[4] * r.rate_coefficient(args...) * Q1(args...)
                        push!(funct_list, f_part)
                        r_counter += 1
                    end

                    # COLLISION INTRINSIC ENERGY GAIN/LOSS
                    r_id = r.reaction_id
                    Er = r.E_threshold             # Collision threshold energy (includes elastic scattering)
                    f_coll(args...) = - Er(args...) * r.rate_coefficient(args...) * args[3] * args[4]
                    push!(funct_list, f_coll)
                    r_counter += 1
                    
                end
            end
            # All the terms for species s are collected.
            # Now, add them to the equation list as a single function 
            f_gainloss(args...) = sum(funct_list[i](args...) for i=1:r_counter)/Q0(args...)
        else
            # If species has no temperature equation then the function returns zero
            f_gainloss(args...) = 0
        end

        # WALL LOSS rates
        # Second, find the WALL LOSS for this specific species
        if (s.has_temp_wall_loss)
            f_wallloss(args...) = s.temp_wall_loss(args...)
        else
            f_wallloss(args...) = 0
        end

        # HEATING MECHANISM
        # Second, find the WALL LOSS for this specific species
        if (s.has_heating_mechanism)
            f_heat(args...) = s.input_power_funct(args...)
        else
            f_heat(args...) = 0
        end

        # This is f in the temperature equation, i.e. in d(T)/dt = f(Te, Tg, n1, n2)
        f_temp(args...) = f_gainloss(args...) + f_wallloss(args...) +
            f_wallloss(args...) + f_heat(args...) 

        # Add f_dens to the list of equations
        push!(temp_eq_list, f_temp)

    end
    return temp_eq_list
end

end