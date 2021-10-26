using InputData: SetupInputData
using InputData: setup_pre_run
using SharedData: c_io_error

function prerun_GM()

    # Reads data from the input.deck file and generated the ReactionSet module
    errcode = SetupInputData("input.deck", setup_pre_run)
    if (errcode == c_io_error)
        return errcode
    end

end