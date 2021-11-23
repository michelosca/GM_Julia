module GM_prerun

using SharedData: c_io_error, Species, Reaction, System, Output, SpeciesID
using InputData_PreRun: SetupInputData!

function prerun_GM(input_file)

    # Reads data from the input.deck file
    species_list = Species[]
    reaction_list = Reaction[]
    system = System()
    speciesID = SpeciesID()
    errcode = SetupInputData!(input_file, species_list, reaction_list, system,
        speciesID)
    if (errcode == c_io_error)
        return errcode
    end

end

end