using SharedData: c_io_error, Species
using InputData: SetupInputData
using InputData: setup_main_run
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using GenerateODEs: GenerateDensRateFunctionList, GenerateTempRateFunctionList
using TestFunctions: TestDensEquations, TestTempEquations, TestRateCoefficients

function run_GM()
    # Reads data from the input.deck file
    errcode, system_list, species_list, reaction_list = SetupInputData("input.deck", setup_main_run)
    if (errcode == c_io_error)
        return errcode
    end

    # Print system, species and reaction lists to terminal
    PrintSystemList(system_list)
    PrintSpeciesList(species_list)
    PrintReactionList(reaction_list)

    dens_funct_list = GenerateDensRateFunctionList(system_list[1], species_list, reaction_list)
    temp_funct_list = GenerateTempRateFunctionList(system_list[1], species_list, reaction_list)

    # Initial conditions
    dens, temp = SetupInitialConditions(species_list)
    print("Initial density values", dens,"\n")
    print("Initial temperature values", temp,"\n")

    for r in reaction_list
        test_val = r.rate_coefficient(temp)
        print("Test val: ", test_val,"\n")
    end

    #TestRateCoefficients(temp, reaction_list)
    TestDensEquations(dens, temp, dens_funct_list, system_list[1])
    TestTempEquations(dens, temp, temp_funct_list, system_list[1])
end

function SetupInitialConditions(species_list::Vector{Species})
    dens = Float64[]
    temp = Float64[]
    for s in species_list
        push!(dens, s.dens0)
        push!(temp, s.temp0)
    end
    return dens, temp
end