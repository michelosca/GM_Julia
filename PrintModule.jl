module PrintModule

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: e, K_to_eV
using SharedData: p_icp_id, p_ccp_id
using InputBlock_Reactions: r_elastic, r_excitat, r_recombi, r_ionizat
using InputBlock_Reactions: r_wall_loss, r_energy_sink, r_cx
using Printf

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - PrintSpeciesList
# - PrintReactionList
# - PrintSystemList


function PrintSpeciesList(species_list::Vector{Species}, sID::SpeciesID)

    @printf("Loaded species\n")
    @printf("%15s %15s %15s %15s %10s %10s %10s %10s %15s %15s\n","Name",
        "Species", "Mass [kg]", "Charge [C]","has n-eq",
        "has T-eq", "has WL", "has P-input","Dens0 [m^-3]","Temp0 [eV]")
    for s in species_list

        # s.neutral_id -> species name
        if (s.species_id == sID.Ar)
            sn_name = "Ar"
        elseif (s.species_id == sID.electron)
            sn_name = "electrons"
        elseif (s.species_id == sID.O2)
            sn_name = "O2"
        elseif (s.species_id == sID.O)
            sn_name = "O"
        else
            sn_name = "None"
        end

        # Density equations
        if (s.has_dens_eq)
            has_dens_eq = "Yes"
        else
            has_dens_eq = "No"
        end

        # Temperature equations
        if (s.has_temp_eq)
            has_temp_eq = "Yes"
        else
            has_temp_eq = "No"
        end

        # Wall losses
        if (s.has_wall_loss)
            has_wall_loss = "Yes"
        else
            has_wall_loss = "No"
        end

        # Power input
        if (s.has_heating_mechanism)
            has_heating = "Yes"
        else
            has_heating = "No"
        end

        @printf("%15s %15s %15.5e %15.5e %10s %10s %10s %10s %15g %15g\n",
            s.name, sn_name, s.mass, s.charge, 
            has_dens_eq, has_temp_eq, has_wall_loss, has_heating,
            s.dens, s.temp*K_to_eV)
    end
    print("\n")
end



function PrintReactionList(reaction_list::Vector{Reaction},
    species_list::Vector{Species}, sID::SpeciesID)

    @printf("Loaded reactions\n")
    @printf("%17s %30s %15s %15s %30s %15s\n",
        "ID - description", "Reaction", "E-threshold [eV]",
        "Neutral","Species balance", "Reactants")

    for r in reaction_list

        r_name = string(r.id)
        # r.case -> reaction name
        if (r.case == r_elastic)
            r_name = string(r_name, "-Elastic")
        elseif (r.case == r_ionizat)
            r_name = string(r_name, "-Ionization")
        elseif (r.case == r_recombi)
            r_name = string(r_name, "-Recombinat.")
        elseif (r.case == r_excitat)
            r_name = string(r_name, "-Excitation")
        elseif (r.case == r_wall_loss)
            r_name = string(r_name, "-Wall react.")
        elseif (r.case == r_energy_sink)
            r_name = string(r_name, "-Energy sink")
        elseif (r.case == r_cx)
            r_name = string(r_name, "-CX")
        end

        # reaction description
        reactants = String[]
        products = String[]
        for i in 1:length(r.involved_species) 
            s = r.involved_species[i]
            b = r.species_balance[i]

            s_name = species_list[s].name
            if b > 0
                for i in 1:b 
                    push!(products, s_name)
                end
            elseif b < 0
                for i in 1:-b
                    push!(reactants, s_name)
                end
            else
                push!(reactants, s_name)
                push!(products, s_name)
            end
        end

        # Reactant string
        react = reactants[1]
        n = length(reactants)
        for i in 2:n
            react = string(react, " + ", reactants[i])
        end
        
        # Product string
        prod = products[1]
        n = length(products)
        for i in 2:n
            prod = string(prod, " + ", products[i])
        end

        # Reaction string
        reaction_str = string(react," -> ",prod)

        # Threshold energy
        E_eV = r.E_threshold / e

        # Neutral species
        r_neutral = ""
        for id in r.neutral_species_id
            r_neutral = string(r_neutral, species_list[id].name,", ")
        end
        r_neutral = chop(r_neutral, tail=2)

        # Involved species and its balance
        r_involved = ""
        i = 1
        for id in r.involved_species
            r_involved = string(r_involved, species_list[id].name,":",r.species_balance[i],", ")
            i += 1
        end
        r_involved= chop(r_involved, tail=2)

        # Reacting species
        r_reactants= ""
        for id in r.reactant_species
            r_reactants = string(r_reactants, species_list[id].name,", ")
        end
        r_reactants = chop(r_reactants, tail=2)

        @printf("%17s %30s %15.2f %15s %30s %15s\n", r_name, reaction_str,
            E_eV, r_neutral, r_involved, r_reactants)
    end
    print("\n")
end


function PrintSystemList(s::System)

    @printf("System parameters\n")
    @printf(" - Area:               %15g m^2\n", s.A)
    @printf(" - Volume:             %15g m^3\n", s.V)
    @printf(" - Length:             %15g m\n", s.l)
    if (s.power_input_method == p_icp_id)
        power_str = "ICP"
    elseif (s.power_input_method == p_ccp_id)
        power_str = "CCP"
    else
        power_str = "Not defined"
    end
    print(" - Input power method: ",power_str,"\n")
    @printf(" - Electrode area:     %15g m^2\n", s.electrode_area)
    @printf(" - Driving frequency:  %15.2f MHz\n", s.drivf)
    @printf("                       %15.f rad/s\n", s.drivOmega)
    @printf(" - Driving power:     %15.2f W\n", s.drivP)
    @printf(" - Driving voltage:   %15.2f V\n", s.drivV)
    @printf(" - Driving current:   %15.2f A\n", s.drivI)
    @printf(" - Simulation time:   %15g s\n", s.t_end)
    print("\n")
end

end