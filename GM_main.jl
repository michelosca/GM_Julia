module GM_main

using SharedData: c_io_error, Species, Reaction, System, SpeciesID, OutputBlock
using InputData: SetupInputData!
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using OutputModule: GenerateOutputs!
using CSV

function run_GM(input)

    if typeof(input) == String
        # Reads data from the input.deck file
        species_list = Species[]
        reaction_list = Reaction[]
        system = System()
        output_list = OutputBlock[]
        speciesID = SpeciesID()
        errcode = SetupInputData!(input, species_list,
            reaction_list, system, output_list, speciesID)
        if (errcode == c_io_error) return errcode end

    elseif typeof(input) == Tuple{Vector{Species}, Vector{Reaction}, System,
        Vector{OutputBlock}, SpeciesID} 
        # Restart from previous simulation
        species_list = input[1]
        reaction_list = input[2]
        system = input[3]
        output_list = input[4]
        speciesID = input[5]
    else
        print("***ERROR*** Input is not recognized\n")
        return 
    end

    # Print system, species and reaction lists to terminal
    PrintSystemList(system)
    PrintSpeciesList(species_list, speciesID)
    PrintReactionList(reaction_list, species_list, speciesID)

    errcode = @time GenerateOutputs!(species_list, reaction_list, system,
        output_list, speciesID)
    if (errcode == c_io_error) return errcode end

    for output in output_list
        label = ""
        for i in 1:output.n_parameters
            label = string(label,"_vs_",output.name[i])
        end
        T_filename = string("T",label,".csv")
        n_filename = string("n",label,".csv")
        K_filename = string("K",label,".csv")
        CSV.write(T_filename, output.T_data_frame )
        CSV.write(n_filename, output.n_data_frame )
        CSV.write(K_filename, output.K_data_frame )
    end

    return species_list, reaction_list, system, speciesID, output_list
end

end