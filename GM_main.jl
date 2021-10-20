#push!(LOAD_PATH, pwd())
using SharedData
using ReadInput: SetupInputData
using GenerateODEs: GenerateDensRateFunctionList, GenerateTempRateFunctionList
using GenerateODEs: GatherListOfFunctions
using Printf
using WallFluxModels: DensWallFluxFunction

global dens = Float64[]
global temp = Float64[]

function run_GM()
    # Reads data from the input.deck file
    errcode = SetupInputData()
    if (errcode != 0)
        return errcode
    end
    PrintSpeciesList()
    PrintReactionList()

    dens_funct_list = GenerateDensRateFunctionList()
    temp_funct_list = GenerateTempRateFunctionList()

    # Initial conditions
    # Density
    n_e = 1.e14
    n_Ar = 1.e20
    n_ArIon = n_e
    dens = [n_e, n_Ar, n_ArIon]

    # Temp
    T_e = 3 * SharedData.e / SharedData.kb
    T_Ar = 300
    T_ArIon = T_Ar
    temp = [T_e, T_Ar, T_ArIon]
    return 1
end

function PrintSpeciesList()

    @printf("Species\n")
    @printf("%15s %15s %15s %15s %15s %15s %15s %15s %15s\n","Name",
        "Species", "Elastic ID", "Mass [kg]", "Charge [C]","has n-eq",
        "has T-eq", "has WL", "has P-input")
    for s in SharedData.species_list

        # s.id -> species name
        if (s.id == SharedData.s_electron_id)
            s_name = "electrons"
        elseif (s.id == SharedData.s_Ar_id)
            s_name = "Ar"
        elseif (s.id == SharedData.s_ArIon_id)
            s_name = "Ar+"
        else
            s_name = "Not found!"
        end

        # s.neutral_id -> species name
        if (s.species_id == SharedData.s_Ar_id)
            sn_name = "Ar"
        elseif (s.species_id == SharedData.s_electron_id)
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

        @printf("%15s %15s %15i %15.5e %15.5e %15s %15s %15s %15s\n",
            s_name, sn_name, s.r_elastic_id, s.mass, s.charge, 
            has_dens_eq, has_temp_eq, has_wall_loss, has_heating)
    end
end



function PrintReactionList()

    @printf("Reactions\n")
    @printf("%15s %30s %15s %15s\n","Name", "Reaction", "E-threshold [eV]",
        "Neutral species")
    for r in SharedData.reaction_list

        # r.id -> reaction name
        if (r.id == SharedData.r_elastic_id)
            r_name = "Elastic"
        elseif (r.id == SharedData.r_ionizat_id)
            r_name = "Ionization"
        elseif (r.id == SharedData.r_recombi_id)
            r_name = "Recombination" 
        elseif (r.id == SharedData.r_excitat_id)
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

            if (s == SharedData.s_electron_id)
                s_name = "e"
            elseif (s == SharedData.s_Ar_id)
                s_name = "Ar"
            elseif (s == SharedData.s_ArIon_id)
                s_name = "Ar+"
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
        E_eV = r.E_threshold / SharedData.e

        # Neutral species
        if (r.neutral_species_id == SharedData.s_Ar_id)
            r_neutral = "Ar"
        else
            r_neutral = "Not found!"
        end

        @printf("%15s %30s %15.2f %15s\n", r_name, reaction_str, E_eV, r_neutral)
    end
end