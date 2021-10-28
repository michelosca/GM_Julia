module GM_main

using SharedData: c_io_error, Species
using InputData: SetupInputData
using InputData: setup_main_run
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using GenerateODEs: GenerateDensRateFunctionList, GenerateTempRateFunctionList
#using TestFunctions: TestDensEquations, TestTempEquations, TestRateCoefficients

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


    #print("Test Rate Coefficients\n")
    #TestRateCoefficients(temp, reaction_list)
    #print("Test Density equations\n")
    #TestDensEquations(dens, temp, dens_funct_list, system_list[1])
    #print("Test Temp equations\n")
    #TestTempEquations(dens, temp, temp_funct_list, system_list[1])
    
    eq_list = cat(temp_funct_list, dens_funct_list,dims=1)
    init = cat(temp,dens,dims=1)
    return init, eq_list
end

function SetupInitialConditions(species_list::Vector{Species})
    dens = Float64[]
    temp = Float64[]
    id = 1
    for s in species_list
        if (s.has_dens_eq)
            push!(dens, s.dens0)
        end
        if (s.has_temp_eq)
            push!(temp, s.temp0)
        end

        id += 1
    end
    return dens, temp
end

end