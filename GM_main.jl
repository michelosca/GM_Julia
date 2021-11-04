module GM_main

using SharedData: c_io_error, Species
using InputData: SetupInputData
using InputData: setup_main_run
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using GenerateODEs: GenerateDensRateFunctionList, GenerateTempRateFunctionList
#using TestFunctions: TestDensEquations, TestTempEquations, TestRateCoefficients
using SolveSystem: ExecuteProblem

function run_GM(input_file)
    # Reads data from the input.deck file
    errcode, system_list, species_list, reaction_list = SetupInputData(input_file, setup_main_run)
    if (errcode == c_io_error)
        return errcode
    end

    # Print system, species and reaction lists to terminal
    PrintSystemList(system_list)
    PrintSpeciesList(species_list)
    PrintReactionList(reaction_list)

    system = system_list[1]
    dens_funct_list = GenerateDensRateFunctionList(system, species_list, reaction_list)
    temp_funct_list = GenerateTempRateFunctionList(system, species_list, reaction_list)

    # Initial conditions
    dens, temp = SetupInitialConditions(species_list)

    #print("Test Rate Coefficients\n")
    #TestRateCoefficients(temp, reaction_list)
    #print("Test Density equations\n")
    #TestDensEquations(dens, temp, dens_funct_list, system_list[1])
    #print("Test Temp equations\n")
    #TestTempEquations(dens, temp, temp_funct_list, system_list[1])
    
    list_offset = length(temp_funct_list)
    eq_list = cat(temp_funct_list, dens_funct_list,dims=1)
    init = cat(temp,dens,dims=1)
    tspan = (0, system.t_end)
    """
    TAr = species_list[2].temp0
    Te = []
    ne = []
    pL = []
    for dens0 in [10^y for y in range(log10(20), log10(24), length=200)]
        kb = 1.380649e-23
        dens[2] = 9.5 * 10^dens0
        init = cat(temp,dens,dims=1)
        p = dens[2] * kb * TAr
        push!(pL, p * system.l)
        print("P*L = ",p*system.l,"; Executing GM problem...\n")
        sol = ExecuteProblem(eq_list, init, tspan, list_offset)
        push!(Te, sol[1,end]/1.16e4)
        push!(ne, sol[5,end])
    end
    return Te, ne, pL 
    """
    sol = ExecuteProblem(eq_list, init, tspan, list_offset)
    return sol
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

end