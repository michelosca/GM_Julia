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
    # Define species ids
    electron_id = 1
    ar_id = 2
    ar_exc_id = 3
    ar_ion_id = 4

    # Define reaction ids
    r_ela_id = 1
    r_exc_id = 2
    r_ion_id = 3
    
    # Define electron-argon elastic scattering
    f_r1(args...) = 2.336e-14 * args[1]^1.609
    Er1(args...) = 3 * me / mg * kb * (args[1] - args[2])
    r1 = Reaction(r_ela_id, [electron_id,ar_id],[0,0], f_r1, Er1)
    # Define e-impact excitation: e + Ar -> e + Ar*
    f_r2(args...) = 2.48e-14 * args[1]^0.33 * exp(-12.78/args[1])
    Er2(args...) = 12.14 * e
    r2 = Reaction(r_exc_id, [electron_id,ar_id,ar_exc_id],[0,-1,1], f_r2, Er2)
    # Define e-impact ionization: e + Ar -> e + e + Ar+
    f_r3(args...) = 2.34e-14 * args[1]^0.59 * exp(-17.44/args[1])
    Er3(args...) =  15.76 * e
    r3 = Reaction(r_ion_id, [electron_id,ar_id,ar_ion_id],[1,-1,1], f_r3, Er3)

    # Load species
    f_zero(args...) = 0
    #                   species_id ,Ln_eq,  LT_eq, Bwln , fwln,                  BwlT, fwlT
    f_dens_wallloss1(args...) = 1 
    f_temp_wallloss1(args...) = 1
    electrons = Species(electron_id, true,  true , true , f_dens_wallloss1, true , f_temp_wallloss1, true,  f_zero)
    argon =     Species(ar_id      , true,  false, false, f_zero,           false, f_zero,           false, f_zero)
    argon_ex =  Species(ar_exc_id  , true,  false, false, f_zero,           false, f_zero,           false, f_zero)
    f_dens_wallloss2(args...) = 1
    f_temp_wallloss2(args...) = 1
    argon_io =  Species(ar_ion_id  , true,  false, true , f_dens_wallloss2, true , f_temp_wallloss2, false, f_zero)

    species_list = [ electrons, argon, argon_ex, argon_io ]
    reaction_list = [r1,r2,r3]
    print("Species type: ",typeof(species_list), "\nReaction type: ", typeof(reaction_list), "\n")
    dens_funct_list = GenerateODEs.GenerateDensODEs(species_list, reaction_list)
    temp_funct_list = GenerateODEs.GenerateTempODEs(species_list, reaction_list)

    return dens_funct_list
end