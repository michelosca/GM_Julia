include("SharedData.jl")
using .SharedData

include("ReadInput.jl")
using .ReadInput

include("GenerateODEs.jl")
using .GenerateODEs

function run_GM()
    GetInputData("input.deck")
    species_list, reaction_list = InitializeData()

    dens_funct_list = GenerateDensRateFunctionList(species_list, reaction_list)
    temp_funct_list = GenerateTempRateFunctionList(species_list, reaction_list)

    ne = 1.e14
    n0 = 1.e20
    dens = [ne, n0, ne]
    Te_eV = 3
    Tg_eV = 0.01
    temp = [Te_eV, Tg_eV]

    dens_val = GatherListOfFunctions(dens_funct_list[1], dens, temp) 
    print("Dens values: ", dens_val,"\n")
    temp_val = GatherListOfFunctions(temp_funct_list[1], dens, temp) 
    print("Temp values: ", temp_val,"\n")

    #test_electron_dens(args::Vector{Float64}) = args[3]*args[4]*f_r3(args) - args[3]*args[5]*f_r4(args)
    #arguments = [3, 0.01, 1.e14, 1.e20, 1.e14]
    #test_val = test_electron_dens(arguments)
    #print("Test values: ", test_val,"\n")

    return 
end