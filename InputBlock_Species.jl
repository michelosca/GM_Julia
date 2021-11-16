module InputBlock_Species

using SharedData: K_to_eV, e, me, amu 
using SharedData: c_io_error
using SharedData: Species

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
        global input_n_id= 0 
        global input_mass = 0.0
        global input_charge = 0.0
        global input_neq_flag = false
        global input_Teq_flag = false
        global input_wl_flag = false
        global input_P_flag = false
        global input_dens0 = 0.0
        global input_temp0 = 0.0
    end
    return errcode
end


function ReadSpeciesEntry(name, var, read_step)

    errcode = c_io_error 

    if (name=="name")
        if (read_step == 2)
            if (var=="e" || var=="electrons" || var=="electron")
                global input_n_id = s_electron_id
                errcode = 0
            elseif (occursin("Ar",var))
                global input_n_id = s_Ar_id
                errcode = 0
            elseif (occursin("O2",var))
                global input_n_id = s_O2_id
                errcode = 0
            elseif (occursin("O",var))
                global input_n_id = s_O_id
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
            global input_charge = parse(Int64,var) * e
        end
        errcode = 0
    end

    if (name=="mass")
        if (read_step == 2)
            expr = Meta.parse(var)
            global input_mass = eval(expr)
        end
        errcode = 0
    end

    if (name=="solve_dens")
        if (read_step == 2)
            global input_neq_flag = parse(Bool, var) 
        end
        errcode = 0
    end

    if (name=="solve_temp")
        if (read_step == 2)
            global input_Teq_flag = parse(Bool, var)
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
            global input_temp0 = parse(Float64, var)*units
        end
        errcode = 0
    end

    if (name=="density" || name=="dens")
        if (read_step == 2)
            global input_dens0 = parse(Float64, var)
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
        if (input_charge != 0)
            if (input_Teq_flag)
                global input_P_flag = true
            end
            global input_wl_flag = true
        end
        
        if (input_mass == 0)
            print("***ERROR*** Species mass has not been defined")
            return c_io_error
        end

        if (input_temp0 == 0)
            print("***ERROR*** Species temperature has not been defined")
            return c_io_error
        end

        errcode = AddSpeciesToList(input_id, input_mass,
            input_charge, input_neq_flag, input_Teq_flag, input_wl_flag,
            input_P_flag, input_n_id, input_dens0, 
            input_temp0)
    end
    return errcode 
end


function AddSpeciesToList(id::Int64, mass::Float64, charge::Float64,
    neq_flag::Bool, Teq_flag::Bool, wl_flag::Bool, P_flag::Bool,
    n_id::Int64, dens0::Float64, temp0::Float64)
    
    errcode = 0
    try
        # Add react to reaction_list
        species = Species(id, n_id, mass, charge, neq_flag,
            Teq_flag, wl_flag, P_flag, dens0, temp0) 
        push!(species_list, species)
    catch
        print("***ERROR*** While attaching species\n")
        errcode = c_io_error
    end

    return errcode
end

end