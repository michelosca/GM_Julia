module ReadInput

using SharedData: me, mp, kb, e
using SharedData: r_elastic_id, r_ionizat_id, r_recombi_id
using SharedData: Reaction, Species
using SharedData: SetupSystemParameters
using SharedData: AddSpeciesToList, AddReactionToList
using SharedData: SetSpeciesID
using SharedData: s_electron_id, s_Ar_id, s_ArIon_id

using PowerInput: PowerInputFunction

# Constants
const c_io_block_end = 0::Int64
const c_io_block_start = 1::Int64
const c_io_entry = 2::Int64
const c_io_emptyline = 3::Int64

function SetupInputData()

    errcode = GetInputData("input.deck")
    if (errcode != 0)
        print("***ERROR***\nFailed to read the input deck. Abort code.")
        return errcode 
    end
    SetupSystemParameters()

   # # Load species: ELECTRONS
   # id = 1
   # SetSpeciesID("e", id)
   # n_id = s_electron_id
   # mass = me
   # charge = -e
   # neq_flag = true
   # Teq_flag = true
   # wl_flag = true
   # P_flag = true
   # relastic_id = 1 
   # AddSpeciesToList(id, mass, charge, neq_flag, Teq_flag, wl_flag, P_flag,
   #     n_id, relastic_id)

   # # Load species: ARGON neutral 
   # id += 1
   # SetSpeciesID("Ar", id)
   # n_id = s_Ar_id
   # mass = 4 * mp
   # charge = 0.0
   # neq_flag = true
   # Teq_flag = false
   # wl_flag = false
   # P_flag = false
   # relastic_id = 0
   # AddSpeciesToList(id, mass, charge, neq_flag, Teq_flag, wl_flag, P_flag,
   #     n_id, relastic_id)

   # # Load species: ARGON ions 
   # id += 1
   # SetSpeciesID("Ar+", id)
   # n_id = s_Ar_id
   # mass = 4 * mp
   # charge = e 
   # neq_flag = true
   # Teq_flag = false 
   # wl_flag = true 
   # P_flag = false
   # relastic_id = 2
   # AddSpeciesToList(id, mass, charge, neq_flag, Teq_flag, wl_flag,
   #     P_flag, n_id, relastic_id)

    ###########################################################################
    # Define electron-argon elastic scattering
    invol_s = [s_electron_id, s_Ar_id]
    balan_s = [0,0]
    react_s = [s_electron_id, s_Ar_id]
    K_r1(temp::Vector{Float64}) = 2.336e-14 * temp[s_electron_id]^1.609
    E = 0.0
    neutral_id = s_Ar_id
    AddReactionToList(r_elastic_id, invol_s, react_s, balan_s, K_r1, E, neutral_id)

    # Define e-impact ionization: e + Ar -> e + e + Ar+
    invol_s = [s_electron_id, s_Ar_id, s_ArIon_id]
    balan_s = [1,-1,1]
    react_s = [s_electron_id, s_Ar_id]
    K_r3(temp::Vector{Float64}) = 2.34e-14 * temp[s_electron_id]^0.59 *
        exp(-17.44/temp[s_electron_id])
    E = 15.76 * e
    neutral_id = s_Ar_id
    AddReactionToList(r_ionizat_id, invol_s, react_s, balan_s, K_r3, E, neutral_id)
    
    # Define Ar recombination: e + Ar+ -> Ar
    invol_s = [s_electron_id, s_ArIon_id, s_Ar_id]
    balan_s = [-1,-1,1]
    react_s = [s_electron_id, s_ArIon_id]
    K_r4(temp::Vector{Float64}) =  5e-39 * temp[s_electron_id]^4.5
    E = 0.0
    neutral_id = s_Ar_id
    AddReactionToList(r_recombi_id, invol_s, react_s, balan_s, K_r4, E, neutral_id)
    return 0
end

global block_id = 0
global b_system = 1
global b_species = 2
global b_reactions = 3

###############################################################################
# START BLOCKS
function GetInputData(filename)

    # Opens the file given in filename and reads each line
    try
        # First read
        print("First reading of the input deck...\n")
        global read_step = 1 
        global input_id = 0
        open(filename,"r") do f
            line = 0
            while ! eof(f)
                s = readline(f)
                # ReadLine identifies each line on filename
                read_flag = ReadLine(s)
                line += 1
            end
        end
        print("End of input deck reading\n\n")

        # Second read
        print("Second reading of the input deck...\n")
        global read_step = 2 
        global input_id = 0
        open(filename,"r") do f
            line = 0
            while ! eof(f)
                s = readline(f)
                # ReadLine identifies each line on filename
                read_flag = ReadLine(s)
                line += 1
            end
        end
        print("End of input deck reading\n\n")
        return 0
    catch
        print("***ERROR***\nFailed reading the input.deck\n")
        return 1
    end

end

function ReadLine(string)

    # Trimm out commets
    i_comment = findfirst("#", string)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        string = string[begin:i_comment-1]
    end

    # Check wether this is a begin/end:block 
    i_block = findfirst(":", string)
    if !(i_block === nothing)
        i_block = i_block[1]
        block_name = string[i_block+1:end]
        if (occursin("begin", string))
            errcode = StartBlock(block_name)
            if (errcode != 0)
                print("***WARNING***\nSomething went wrong starting the ",
                    block_name, " block\n")
            end
            return c_io_block_start
        elseif (occursin("end", string))
            errcode = EndBlock(block_name)
            if (errcode != 0)
                print("***WARNING***\nSomething went wrong ending the ",
                    block_name, " block\n")
            end
            return c_io_block_end
        end
    else
        # Check whether it is a "name = var" line
        i_eq = findfirst("=", string)
        if !(i_eq === nothing)
            i_eq = i_eq[1]
            name = strip(string[begin:i_eq-1])
            var = strip(string[i_eq+1:end])
            #print("Name: ", name, "  ")
            #print(" Var: ", var , "\n\n")
            ReadInputDeckEntry(name, var, block_id)
            return c_io_entry
        end
    end
    return c_io_emptyline
end

###############################################################################
# BLOCK START/END
function StartBlock(name)

    errcode = 0
    if (occursin("system",name))
        global block_id = b_system
        errcode = StartSystemBlock()
    elseif (occursin("species",name))
        global block_id = b_species
        errcode = StartSpeciesBlock()
    elseif (occursin("reactions",name))
        global block_id = b_reactions
        errcode = StartReactionsBlock()
    end
    return errcode
end

function EndBlock(name)

    errcode = 0
    global block_id = 0
    if (occursin("system",name))
        errcode = EndSystemBlock()
    elseif (occursin("species",name))
        errcode = EndSpeciesBlock()
    elseif (occursin("reactions",name))
        errcode = EndReactionsBlock()
    end
    return errcode
end
#
###############################################################################
###############################################################################
# START BLOCKS
function StartSystemBlock()
    return 0
end


function StartSpeciesBlock()

    global input_id += 1
    if (read_step == 2)
        global input_n_id = 0
        global input_relastic_id = 0 
        global input_mass = 0.0
        global input_charge = 0.0
        global input_neq_flag = false 
        global input_Teq_flag = false 
        global input_wl_flag = false
        global input_P_flag = false
    end
    return 0
end

function StartReactionsBlock()
    return 0
end
#
###############################################################################
###############################################################################
# END BLOCKS
function EndSystemBlock()
    return 0
end

function EndSpeciesBlock()

    if (read_step == 1)
        return 0
    elseif (read_step == 1)
        if (input_charge != 0)
            if (input_Teq_flag)
                global input_P_flag = true
            end
            global input_wl_flag = true
        end

        AddSpeciesToList(input_id, input_mass, input_charge, input_neq_flag,
            input_Teq_flag, input_wl_flag, input_P_flag, input_n_id, input_relastic_id)
        return 0
    end
end

function EndReactionsBlock()
    return 0
end
#
###############################################################################
###############################################################################
# READ INPUT ENTRIES 
function ReadInputDeckEntry(name, var, block_id)

    errcode = 0
    
    if (block_id == b_system)
        errcode = ReadSystemEntry(name, var)
    elseif (block_id == b_species)
        errcode = ReadSpeciesEntry(name, var)
    elseif (block_id == b_reactions)
        errcode = ReadReactionsEntry(name, var)
    end

    if (errcode == 1)
        print("***WARNING***\nEntry in System block has not beeb identified\n")
        print("Entry ", string(name,"=",var) ," does not belong to any block\n")
    end
   
end

function ReadSystemEntry(name, var)
    units_index = findlast("_", var)
    units_index = units_index[1]
    units = var[units_index+1:end]
    var = var[1:units_index-1]
    
    errcode = SharedData.SetSystemParameters(name, var, unis)
    return errcode 
end

function ReadSpeciesEntry(name, var)

    if (name=="name")
        if (read_step == 1)
            errcode = SetSpeciesID(var, input_id)
        else 
            if (var=="e" || var=="electrons" || var=="electron")
                global input_n_id = s_electron_id
            elseif (occursin("Ar",var))
                global input_n_id = s_Ar_id
            else
                print("***WARNING***\n")
                print("Neutral species id has not been found. Make sure you\n",
                    "define first electron and neutral species.")
            end
        end
        return 0
    end

    if (read_step == 2)
        if (name=="charge")
            global input_charge = parse(Int64,var) * e
            return 0
        end

        if (name=="mass")
            expr = Meta.parse(var)
            global input_mass = eval(expr)
            return 0
        end

        if (name=="solve_dens")
            global input_neq_flag = parse(Bool, var) 
            return 0
        end

        if (name=="solve_temp")
            global input_Teq_flag = parse(Bool, var)
            return 0
        end
    end
    return 1
end

function ReadReactionsEntry(name, var)
    return 0 
end

end