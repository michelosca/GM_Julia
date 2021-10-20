module ReadInput

using SharedData: SetSystemParameters
using SharedData: me, mp, kb, e
using SharedData: r_elastic_id, r_ionizat_id, r_recombi_id
using SharedData: Reaction, Species
using SharedData: AddSpeciesToList, AddReactionToList
using SharedData: SetSpeciesID
using SharedData: s_electron_id, s_Ar_id, s_ArIon_id

using PowerInput: PowerInputFunction

# Constants
const c_io_block_end = 0::Int64
const c_io_error = 1::Int64
const c_io_block_start = 4::Int64
const c_io_entry = 2::Int64
const c_io_emptyline = 3::Int64

function SetupInputData()

    errcode = ReadInputData("input.deck")
    if (errcode != 0)
        print("***ERROR*** Failed to read the input deck. Abort code\n")
        return errcode 
    end

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
function ReadInputData(filename)

    errcode = 1

    # Opens the file given in filename and reads each line
    # First read
    print("First reading of the input deck...\n")
    global read_step = 1 
    global input_id = 0
    errcode = ReadFile(filename)
    if (errcode == c_io_error)
        return errcode
    end
    print("End of input deck reading\n\n")

    # Second read
    print("Second reading of the input deck...\n")
    global read_step = 2 
    global input_id = 0
    errcode = ReadFile(filename)
    if (errcode == c_io_error)
        return errcode
    end
    print("End of input deck reading\n\n")
    return errcode
end


function ReadFile(filename::String)
    
    errcode = 0

    open(filename,"r") do f
        line = 1
        while ! eof(f)
            s = readline(f)
            # ReadLine identifies each line on filename
            errcode = ReadLine(s)
            if (errcode == c_io_error)
                print("***ERROR*** Stop reading at file line ", line,"\n")
                return errcode  
            end
            line += 1
        end
    end
    return errcode
end


function ReadLine(string)

    errcode = 1

    # Trimm out commets
    i_comment = findfirst("#", string)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        string = string[begin:i_comment-1]
    end

    # Signal wether this is a begin/end:block 
    i_block = findfirst(":", string)
    # Signla whether it is a "name = var" line
    i_eq = findfirst("=", string)

    # Check line
    if !(i_block === nothing)
        i_block = i_block[1]
        global block_name = string[i_block+1:end]
        if (occursin("begin", string))
            errcode = StartBlock(block_name)
            if (errcode == c_io_error)
                print("***WARNING***\n  - Something went wrong starting the ",
                    block_name, " block\n")
            end
        elseif (occursin("end", string))
            errcode = EndBlock(block_name)
            if (errcode == c_io_error)
                print("***WARNING***\n  - Something went wrong ending the ",
                    block_name, " block\n")
            end
        end
    elseif !(i_eq === nothing)
        i_eq = i_eq[1]
        name = strip(string[begin:i_eq-1])
        var = strip(string[i_eq+1:end])
        errcode = ReadInputDeckEntry(name, var, block_id)
        if (errcode == c_io_error)
            print("***WARNING***\n  - Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", name ," = ",var ,"\n")
        end
    else
        errcode = c_io_emptyline
    end
    return errcode 
end

###############################################################################
# BLOCK START/END
function StartBlock(name)

    errcode = c_io_error
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

    errcode = c_io_error
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
    errcode = 0
    return errcode
end


function StartSpeciesBlock()

    errcode = 0
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
    return errcode
end

function StartReactionsBlock()
    errcode = 0
    return errcode
end
#
###############################################################################
###############################################################################
# END BLOCKS
function EndSystemBlock()
    errcode = 0
    return errcode
end

function EndSpeciesBlock()

    errcode = c_io_error 

    if (read_step == 1)
        errcode = 0
    elseif (read_step == 2)
        if (input_charge != 0)
            if (input_Teq_flag)
                global input_P_flag = true
            end
            global input_wl_flag = true
        end

        errcode = AddSpeciesToList(input_id, input_mass, input_charge, input_neq_flag,
            input_Teq_flag, input_wl_flag, input_P_flag, input_n_id, input_relastic_id)
    end
    return errcode 
end

function EndReactionsBlock()
    errcode = 0
    return errcode
end
#
###############################################################################
###############################################################################
# READ INPUT ENTRIES 
function ReadInputDeckEntry(name, var, block_id)

    errcode = 1 
    
    if (block_id == b_system)
        errcode = ReadSystemEntry(name, var)
    elseif (block_id == b_species)
        errcode = ReadSpeciesEntry(name, var)
    elseif (block_id == b_reactions)
        errcode = ReadReactionsEntry(name, var)
    end

    return errcode
end

function ReadSystemEntry(name, var)

    errcode = 1

    if (read_step == 1)
        units_fact = 1.0
        units_index = findlast("_", var)
        if (units_index === nothing)
            units = ""
        else
            units_index = units_index[1]
            units_str = lowcase(var[units_index+1:end])
            if (units_str=="kw" || units_str=="khz")
                units_fact = 1.e3
                var = var[1:units_index-1]
            elseif (units_str=="mhz" || units_str=="mw")
                units_fact = 1.e6
                var = var[1:units_index-1]
            elseif (units_str=="ghz")
                units_fact = 1.e9
                var = var[1:units_index-1]
            end
        end
        
        errcode = SetSystemParameters(name, var, units_fact)
    else
        errcode = 0
    end
    return errcode 
end

function ReadSpeciesEntry(name, var)

    errcode = c_io_error 

    if (name=="name")
        if (read_step == 1)
            errcode = SetSpeciesID(var, input_id)
        else 
            if (var=="e" || var=="electrons" || var=="electron")
                global input_n_id = s_electron_id
                errcode = 0
            elseif (occursin("Ar",var))
                global input_n_id = s_Ar_id
                errcode = 0
            else
                print("***WARNING***\n")
                print("Neutral species id has not been found. Make sure you\n",
                    "define first electron and neutral species.")
                errcode = c_io_error 
            end
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
    return errcode 
end

function ReadReactionsEntry(name, var)
    errcode = 0
    return errcode
end

end