module InputData

using SharedData: c_io_error

using InputBlock_Species: StartFile_Species
using InputBlock_Species: StartSpeciesBlock, EndSpeciesBlock, ReadSpeciesEntry
using InputBlock_Species: species_list

using InputBlock_Reactions: StartFile_Reactions
using InputBlock_Reactions: StartReactionsBlock, EndReactionsBlock, ReadReactionsEntry
using InputBlock_Reactions: reaction_list

using InputBlock_System: StartFile_System
using InputBlock_System: StartSystemBlock, EndSystemBlock, ReadSystemEntry
using InputBlock_System: system_list

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
const setup_pre_run = 0
const setup_main_run = 1

# INPUT BLOCK IDs
global block_id = 0
const b_system = 1
const b_species = 2
const b_reactions = 3

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - SetupInputData
#   - ReadInputData
#     - StartFile
#     - ReadFile
#       - ReadLine
#         - StartBlock
#           - StartSystemBlock
#           - StartSpeciesBlock
#           - StartReactionsBlock
#         - ReadInputDeckEntry
#           - ReadSystemEntry
#           - ReadSpeciesEntry
#           - ReadReactionsEntry
#         - EndBlock
#           - EndSystemBlock
#           - EndtSpeciesBlock
#           - EndtReactionsBlock

function SetupInputData(filename, setup_flag)

    errcode = c_io_error

    errcode = ReadInputData(filename, setup_flag)
    if (errcode == c_io_error)
        print("***ERROR*** Failed to read the input deck. Abort code\n")
    end

    if (setup_flag == setup_pre_run)
        return errcode
    elseif (setup_flag == setup_main_run)
        return errcode, system_list, species_list, reaction_list 
    end
end


function ReadInputData(filename, setup_flag)

    errcode = c_io_error 

    if (setup_flag == setup_pre_run)
        read_step = 0 
        print("Reading $read_step of the input deck...\n")
        errcode = StartFile(read_step)
        if (errcode == c_io_error) return errcode end
        errcode = ReadFile(filename, read_step)
        if (errcode == c_io_error) return errcode end
        print("End of input deck reading\n\n")

    elseif (setup_flag == setup_main_run)
        # Opens the file given in filename and reads each line
        for read_step in 1:2
            print("Reading $read_step of the input deck...\n")
            errcode = StartFile(read_step)
            if (errcode == c_io_error) return errcode end
            errcode = ReadFile(filename, read_step)
            if (errcode == c_io_error) return errcode end
            print("End of input deck reading\n\n")
        end
    else
        print("***ERROR*** Setup flag is not recognized\n")
    end
    return errcode
end


function ReadFile(filename::String, read_step::Int64)
    
    errcode = 0

    open(filename,"r") do f
        line = 1
        while ! eof(f)
            s = readline(f)
            # ReadLine identifies each line on filename
            errcode = ReadLine(s, read_step)
            if (errcode == c_io_error)
                print("***ERROR*** Stop reading at file line ", line,"\n")
                return errcode  
            end
            line += 1
        end
    end
    return errcode
end


function ReadLine(string, read_step)

    errcode = c_io_error 

    # Trimm out commets
    i_comment = findfirst("#", string)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        string = string[begin:i_comment-1]
    end
    string = strip(string)

    # Signal wether this is a begin/end:block 
    i_block = findfirst(":", string)
    # Signla whether it is a "name = var" line
    i_eq = findfirst("=", string)
    # Signla whether it is a reaction line
    i_react = findfirst(";", string)

    # Check line
    if !(i_block === nothing)
        i_block = i_block[1]
        global block_name = string[i_block+1:end]
        if (occursin("begin", string))
            errcode = StartBlock(block_name, read_step)
            if (errcode == c_io_error)
                print("***ERROR*** Something went wrong starting the ",
                    block_name, " block\n")
                return c_io_error
            end
        elseif (occursin("end", string))
            errcode = EndBlock(block_name, read_step)
            if (errcode == c_io_error)
                print("***ERROR*** Something went wrong ending the ",
                    block_name, " block\n")
                return c_io_error
            end
        end
    elseif !(i_eq === nothing)
        i_eq = i_eq[1]
        name = strip(string[begin:i_eq-1])
        var = strip(string[i_eq+1:end])
        errcode = ReadInputDeckEntry(name, var, block_id, read_step)
        if (errcode == c_io_error)
            print("***WARNING*** Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", name ," = ",var ,"\n")
        end
    elseif !(i_react === nothing)
        errcode = ReadInputDeckEntry("",string, block_id, read_step) 
        if (errcode == c_io_error)
            print("***WARNING*** Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", string ,"\n")
        end
    else
        errcode = 0
    end
    return errcode 
end


function StartFile(read_step)

    errcode = StartFile_Species(read_step) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input species block")
        return errcode
    end

    errcode = StartFile_Reactions(read_step) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input reaction block")
        return errcode
    end
    
    errcode = StartFile_System(read_step) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input system block")
        return errcode
    end
    
    return errcode
end


function StartBlock(name, read_step)

    errcode = c_io_error
    if (occursin("system",name))
        global block_id = b_system
        errcode = StartSystemBlock(read_step)
    elseif (occursin("species",name))
        global block_id = b_species
        errcode = StartSpeciesBlock(read_step)
    elseif (occursin("reactions",name))
        global block_id = b_reactions
        errcode = StartReactionsBlock(read_step)
    end
    return errcode
end

function EndBlock(name, read_step)

    errcode = c_io_error
    global block_id = 0
    if (occursin("system",name))
        errcode = EndSystemBlock(read_step)
    elseif (occursin("species",name))
        errcode = EndSpeciesBlock(read_step)
    elseif (occursin("reactions",name))
        errcode = EndReactionsBlock(read_step)
    end
    return errcode
end


function ReadInputDeckEntry(name, var, block_id, read_step)

    errcode = 1 
    
    if (block_id == b_system)
        errcode = ReadSystemEntry(name, var, read_step)
    elseif (block_id == b_species)
        errcode = ReadSpeciesEntry(name, var, read_step)
    elseif (block_id == b_reactions)
        errcode = ReadReactionsEntry(name, var, read_step)
    end

    return errcode
end

end