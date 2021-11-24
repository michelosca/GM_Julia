module GM_main

using SharedData: c_io_error, Species, Reaction, System, Output, SpeciesID
using InputData: SetupInputData!
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using OutputModule: GenerateOutputs

function run_GM(input, output_flag_list)

    if typeof(input) == String
        # Reads data from the input.deck file
        species_list = Species[]
        reaction_list = Reaction[]
        system = System()
        speciesID = SpeciesID()
        errcode = SetupInputData!(input, species_list,
            reaction_list, system, speciesID)
        if (errcode == c_io_error)
            return
        end
    elseif typeof(input) == Tuple{Vector{Species}, Vector{Reaction}, System,
        SpeciesID} 
        # Restart from previous simulation
        species_list = input[1]
        reaction_list = input[2]
        system = input[3]
        speciesID = input[4]
    else
        print("***ERROR*** Input is not recognized\n")
        return 
    end

    # Print system, species and reaction lists to terminal
    PrintSystemList(system)
    PrintSpeciesList(species_list, speciesID)
    PrintReactionList(reaction_list, species_list, speciesID)

    # Execute problem(s)
    output_list = Output[]
    output = Output(output_flag_list)
    push!(output_list, output)

    output_list = @time GenerateOutputs(species_list, reaction_list, system,
        speciesID, output_list)

    return species_list, reaction_list, system, speciesID, output_list
end

end