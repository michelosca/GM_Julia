include("SharedData.jl")
using .SharedData

include("ReadInput.jl")
using .ReadInput

include("GenerateODEs.jl")
using .GenerateODEs

# argms: function arguments
# [1] electron temperature in eV
# [2] neutral gas temperature in eV
# [3] species 1 number density
# [4] species 2 number density


function run_GM()
    ReadInput.GetInputData("input.deck")

    #dens_funct_list = DensityODEs.Initialize()
    
    ar_id = 2
    ar_ion_id = 3
    ar_exc_id = 4

    # Define reaction ids
    r_ela_id = 1
    r_exc_id = 2
    r_ion_id = 3
    r_rec_id = 4

    # Define electron-argon elastic scattering
    f_r1(temp::Vector{Float64}) = 2.336e-14 * temp[c_electron_id]^1.609
    Er1(temp::Vector{Float64}) = 3 * me / mg * kb * (temp[c_electron_id] - temp[c_neutral_id])
    r1 = Reaction(r_ela_id, [c_electron_id,ar_id],[c_electron_id, ar_id],[0,0], f_r1, Er1)
    # Define e-impact ionization: e + Ar -> e + e + Ar+
    f_r3(temp::Vector{Float64}) = 2.34e-14 * temp[c_electron_id]^0.59 * exp(-17.44/temp[c_electron_id])
    Er3(temp::Vector{Float64}) =  15.76 * e
    r3 = Reaction(r_ion_id, [c_electron_id,ar_id,ar_ion_id],[c_electron_id, ar_id], [1,-1,1], f_r3, Er3)
    # Define Ar recombination: e + Ar+ -> Ar
    f_r4(temp::Vector{Float64}) =  5e-39 * temp[c_electron_id]^4.5
    Er4(temp::Vector{Float64}) =  0
    r4 = Reaction(r_rec_id, [c_electron_id,ar_ion_id, ar_id], [c_electron_id, ar_ion_id], [-1,-1,1], f_r4, Er4)

    # Load species
    f_zero(temp::Vector{Float64}) = 0
    #                   species_id ,Ln_eq,  LT_eq, Bwln , fwln,                  BwlT, fwlT
    f_dens_wallloss1(args::Vector{Float64}) = 1 
    f_temp_wallloss1(args::Vector{Float64}) = 1
    electrons = Species(c_electron_id, true,  true , false , f_dens_wallloss1, false , f_temp_wallloss1, false,  f_zero)
    argon =     Species(ar_id      , true,  false, false, f_zero,           false, f_zero,           false, f_zero)
    argon_ex =  Species(ar_exc_id  , true,  false, false, f_zero,           false, f_zero,           false, f_zero)
    f_dens_wallloss2(args::Vector{Float64}) = 1
    f_temp_wallloss2(args::Vector{Float64}) = 1
    argon_io =  Species(ar_ion_id  , true,  false, true , f_dens_wallloss2, true , f_temp_wallloss2, false, f_zero)

    species_list = [ electrons, argon, argon_io ]
    reaction_list = [r1,r3,r4]
    print("Species type: ",typeof(species_list), "\nReaction type: ", typeof(reaction_list), "\n")
    dens_funct_list = GenerateODEs.GenerateDensRateFunctionList(species_list, reaction_list)
    temp_funct_list = GenerateODEs.GenerateTempRateFunctionList(species_list, reaction_list)

    # args
    # 1 electron temp
    # 2 argon temp
    # 3 electron dens
    # 4 argon dens
    # 5 argon ions dens
    test_electron_dens(args::Vector{Float64}) = args[3]*args[4]*f_r3(args) - args[3]*args[5]*f_r4(args)

    arguments = [3, 0.01, 1.e14, 1.e20, 1.e14]
    ne = 1.e14
    n0 = 1.e20
    dens = [ne, n0, ne]
    Te_eV = 3
    Tg_eV = 0.01
    temp = [Te_eV, Tg_eV]
    code_val = GenerateODEs.gather_list_of_functions(dens_funct_list[1], dens, temp) 
    test_val = test_electron_dens(arguments)
    print("Code values: ", code_val,"\n")
    print("Test values: ", test_val,"\n")

    return 
end