module InputBlock_Species

using SharedData: K_to_eV, e, me, amu 
using SharedData: c_io_error
using SharedData: Species, Reaction

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
global species_list = Species[]

# SPECIES IDs 
# In case more species need to be defined, they need to be added here
global s_electron_id = 0

global s_Ar_id = 0
global s_O_id = 0
global s_O2_id = 0

global s_ArIon_id = 0
global s_ArExc_id = 0

global s_OnIon_id = 0
global s_OIon_id = 0
global s_O2Ion_id = 0
global s_O3p_id = 0
global s_O1d_id = 0
global s_O2a1Ag_id = 0

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - StartFile_Species
# - StartSpeciesBlock
# - ReadSpeciesEntry
#   - SetSpeciesID
# - EndSpeciesBlock
#   - AddSpeciesToList


function StartFile_Species(read_step) 

    errcode = 0

    global input_id = 0 
    if (read_step == 2)
        global species_list = Species[]
    end


    return errcode
end

function StartSpeciesBlock(read_step)

    errcode = 0
    global input_id += 1
    if (read_step == 2)
        global current_species = Species()
        current_species.id = input_id
        current_species.species_id = 0
        current_species.mass = 0.0
        current_species.charge = 0.0
        current_species.has_dens_eq = false
        current_species.has_temp_eq = false
        current_species.has_wall_loss = false
        current_species.has_heating_mechanism = false
        current_species.dens = 0.0
        current_species.temp = 0.0
        current_species.reaction_list = Reaction[]
        current_species.mfp = 1.e100
        current_species.v_thermal = 0.0
        current_species.v_Bohm = 0.0
        current_species.D = 0.0
        current_species.Lambda = 0.0
        current_species.h_R = 0.0
        current_species.h_L = 0.0
        current_species.gamma = 0.0
        current_species.n_sheath = 0.0
        current_species.flux = 0.0
    end
    return errcode
end


function ReadSpeciesEntry(name, var, read_step)

    errcode = c_io_error 

    if (name=="name")
        if (read_step == 2)
            if (var=="e" || var=="electrons" || var=="electron")
                current_species.species_id = s_electron_id
                errcode = 0
            elseif (occursin("Ar",var))
                current_species.species_id = s_Ar_id
                errcode = 0
            elseif (occursin("O2",var))
                current_species.species_id = s_O2_id
                errcode = 0
            elseif (occursin("O",var))
                current_species.species_id = s_O_id
                errcode = 0
            else
                print("***ERROR*** Neutral species id has not been found\n")
                errcode = c_io_error 
            end
        else 
            # This is set in the pre-run(read_step==0) and in main-run(read_step===1)
            errcode = SetSpeciesID(var, input_id)
        end
    end

    if (name=="charge")
        if (read_step == 2)
            current_species.charge = parse(Int64,var) * e
        end
        errcode = 0
    end

    if (name=="mass")
        if (read_step == 2)
            expr = Meta.parse(var)
            current_species.mass = eval(expr)
        end
        errcode = 0
    end

    if (name=="solve_dens")
        if (read_step == 2)
            current_species.has_dens_eq = parse(Bool, var) 
        end
        errcode = 0
    end

    if (name=="solve_temp")
        if (read_step == 2)
            current_species.has_temp_eq = parse(Bool, var)
        end
        errcode = 0
    end

    if (name=="T" || name=="T_eV")
        if (read_step == 2)
            if (name=="T_eV")
                units = 1.0/K_to_eV
            else
                units = 1.0
            end
            current_species.temp = parse(Float64, var)*units
        end
        errcode = 0
    end

    if (name=="density" || name=="dens")
        if (read_step == 2)
            current_species.dens = parse(Float64, var)
        end
        errcode = 0
    end
    return errcode 
end


function SetSpeciesID(species_name, id)

    errcode = 0

    if ("e" == species_name || "electrons" == species_name)
        global s_electron_id = id
    elseif ("Ar" == species_name)
        global s_Ar_id = id
    elseif ("O" == species_name)
        global s_O_id = id
    elseif ("O2" == species_name)
        global s_O2_id = id
    elseif ("Ar+" == species_name)
        global s_ArIon_id = id
    elseif ("O+" == species_name)
        global s_OIon_id = id
    elseif ("O-" == species_name)
        global s_OnIon_id = id
    elseif ("O2+" == species_name)
        global s_O2Ion_id = id
    elseif ("Ar*" == species_name || "Ar excited" == species_name)
        global s_ArExc_id = id
    elseif ("O(3p)" == species_name)
        global s_O3p_id = id
    elseif ("O(1d)" == species_name)
        global s_O1d_id = id
    elseif ("O2(a1Ag)" == species_name)
        global s_O2a1Ag_id = id
    else
        errcode = 1
    end
    return errcode 
end


function EndSpeciesBlock(read_step)

    errcode = 0 

    if (read_step == 1)
        errcode = 0
    elseif (read_step == 2)
        if (current_species.charge != 0)
            if (current_species.has_temp_eq)
                current_species.has_heating_mechanism = true
            end
            current_species.has_wall_loss = true
        end
        
        if (current_species.mass == 0)
            print("***ERROR*** Species mass has not been defined")
            return c_io_error
        end

        if (current_species.temp == 0)
            print("***ERROR*** Species temperature has not been defined")
            return c_io_error
        end

        push!(species_list, current_species)
    end
    return errcode 
end

end