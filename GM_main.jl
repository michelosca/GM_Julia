module GM_main

using SharedData: c_io_error, Species, Reaction, System, Output, SpeciesID
using InputData: SetupInputData!, setup_main_run
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using OutputModule: GenerateOutputs

function run_GM(input_file, output_flag_list)
    # Reads data from the input.deck file
    species_list = Species[]
    reaction_list = Reaction[]
    system = System()
    speciesID = SpeciesID()
    errcode = SetupInputData!(input_file, setup_main_run, species_list,
        reaction_list, system, speciesID)
    if (errcode == c_io_error)
        return errcode
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

    return output_list
end

end