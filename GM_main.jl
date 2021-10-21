#push!(LOAD_PATH, pwd())
using SharedData
using PrintModule
using ReadInput: SetupInputData
using GenerateODEs: GenerateDensRateFunctionList, GenerateTempRateFunctionList

global dens = Float64[]
global temp = Float64[]

function run_GM()
    # Reads data from the input.deck file
    errcode = SetupInputData()
    if (errcode != 0)
        return errcode
    end

    # Print species and reaction lists to terminal
    PrintModule.PrintSpeciesList()
    PrintModule.PrintReactionList()

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

    return 
end