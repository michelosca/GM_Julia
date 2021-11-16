module GM_main

using SharedData: c_io_error, Output
using InputData: SetupInputData
using InputData: setup_main_run
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using GenerateODEs: GenerateDensRateFunctionList, GenerateTempRateFunctionList
using OutputModule: GenerateOutputs

function run_GM(input_file, output_flag_list)
    # Reads data from the input.deck file
    errcode, system_list, species_list, reaction_list =
        SetupInputData(input_file, setup_main_run)
    if (errcode == c_io_error)
        return errcode
    end

    # Print system, species and reaction lists to terminal
    PrintSystemList(system_list)
    PrintSpeciesList(species_list)
    PrintReactionList(reaction_list)

    # Generate list with equation terms
    dens_funct_list = GenerateDensRateFunctionList(species_list, reaction_list)
    temp_funct_list = GenerateTempRateFunctionList(species_list, reaction_list)

    # Execute problem(s)
    output_list = Output[]
    output = Output(output_flag_list)
    push!(output_list, output)
    output_list = @time GenerateOutputs(dens_funct_list, temp_funct_list, species_list,
        reaction_list, system_list, output_list)

    return output_list
end

end