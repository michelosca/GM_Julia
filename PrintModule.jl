module PrintModule

using SharedData: System, Species, Reaction
using SharedData: e, K_to_eV
using SharedData: p_icp_id, p_ccp_id
using InputBlock_Species: s_electron_id, s_Ar_id, s_ArIon_id, s_ArExc_id
using InputBlock_Reactions: r_elastic_id, r_excitat_id, r_recombi_id, r_ionizat_id
using Printf

function PrintSpeciesList(species_list::Vector{Species})

    @printf("Loaded species\n")
    @printf("%15s %15s %15s %15s %15s %10s %10s %10s %10s %15s %15s\n","Name",
        "Species", "Elastic ID", "Mass [kg]", "Charge [C]","has n-eq",
        "has T-eq", "has WL", "has P-input","Dens0 [m^-3]","Temp0 [eV]")
    for s in species_list
        # s.id -> species name
        if (s.id == s_electron_id)
            s_name = "electrons"
        elseif (s.id == s_Ar_id)
            s_name = "Ar"
        elseif (s.id == s_ArIon_id)
            s_name = "Ar+"
        elseif (s.id == s_ArExc_id)
            s_name = "Ar*"
        else
            s_name = "Not found!"
        end

        # s.neutral_id -> species name
        if (s.species_id == s_Ar_id)
            sn_name = "Ar"
        elseif (s.species_id == s_electron_id)
            sn_name = "electrons"
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

        @printf("%15s %15s %15i %15.5e %15.5e %10s %10s %10s %10s %15g %15g\n",
            s_name, sn_name, s.r_elastic_id, s.mass, s.charge, 
            has_dens_eq, has_temp_eq, has_wall_loss, has_heating,
            s.dens0, s.temp0*K_to_eV)
    end
    print("\n")
end



function PrintReactionList(reaction_list::Vector{Reaction})

    @printf("Loaded reactions\n")
    @printf("%15s %30s %15s %8s %10s %10s %10s %10s\n",
        "Name", "Reaction", "E-threshold [eV]",
        "Neutral","Involved sp.", "Sp. balance", "Reactant sp.",
        "Rate coeff.")
    for r in reaction_list
        # r.id -> reaction name
        if (r.id == r_elastic_id)
            r_name = "Elastic"
        elseif (r.id == r_ionizat_id)
            r_name = "Ionization"
        elseif (r.id == r_recombi_id)
            r_name = "Recombination" 
        elseif (r.id == r_excitat_id)
            r_name = "Excitation" 
        else
            r_name = "Not found!"
        end

        # reaction description
        reactants = String[]
        products = String[]
        for i in 1:length(r.involved_species) 
            s = r.involved_species[i]
            b = r.species_balance[i]

            if (s == s_electron_id)
                s_name = "e"
            elseif (s == s_Ar_id)
                s_name = "Ar"
            elseif (s == s_ArIon_id)
                s_name = "Ar+"
            elseif (s == s_ArExc_id)
                s_name = "Ar*"
            else
                s_name = "None"
            end
            if b > 0
                push!(products, s_name)
            elseif b < 0
                push!(reactants, s_name)
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
        if (r.neutral_species_id == s_Ar_id)
            r_neutral = "Ar"
        else
            r_neutral = "Not found!"
        end

        @printf("%15s %30s %15.2f %8s %10s %10s %10s %s\n", r_name, reaction_str,
            E_eV, r_neutral,r.involved_species, r.species_balance, r.reactant_species,
            r.rate_coefficient)
    end
    print("\n")
end


function PrintSystemList(system_list::Vector{System})

    @printf("System parameters\n")
    i = 0
    for s in system_list
        i += 1
        print("System $i\n")
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
        @printf(" - Driving frequency:  %15.2f MHz\n", s.drivf)
        @printf("                       %15.f rad/s\n", s.drivOmega)
        @printf(" - Driving power:     %15.2f W\n", s.drivP)
        @printf(" - Driving voltage:   %15.2f V\n", s.drivV)
        @printf(" - Driving current:   %15.2f A\n", s.drivI)
    end
    print("\n")
end

end