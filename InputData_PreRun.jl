module InputData_PreRun

using SharedData: c_io_error
using SharedData: p_icp_id, p_ccp_id 
using SharedData: Species, Reaction, System, SpeciesID

using InputBlock_Species: StartFile_Species!
using InputBlock_Species: StartSpeciesBlock!, EndSpeciesBlock!, ReadSpeciesEntry!

using InputBlock_Reactions_PreRun: StartFile_Reactions!
using InputBlock_Reactions_PreRun: StartReactionsBlock!, EndReactionsBlock!, ReadReactionsEntry!
using InputBlock_Reactions_PreRun: r_elastic, r_wall_loss 

using InputBlock_System: StartFile_System!
using InputBlock_System: StartSystemBlock!, EndSystemBlock!, ReadSystemEntry!

using PlasmaParameters: GetLambda, GetGamma
using EvaluateExpressions: ReplaceSystemSymbolS!

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
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

function SetupInputData!(filename::String, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, speciesID::SpeciesID)

    errcode = c_io_error

    errcode = ReadInputData!(filename, species_list, reaction_list, system,
        speciesID)
    if (errcode == c_io_error)
        print("***ERROR*** Failed to read the input deck. Abort code\n")
    end

    return errcode
end


function ReadInputData!(filename::String, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, speciesID::SpeciesID)

    errcode = c_io_error 

    print("Reading the input deck...\n")
    read_step = 0
    errcode = StartFile!(read_step, species_list, reaction_list, system,
        speciesID)
    if (errcode == c_io_error) return errcode end
    errcode = ReadFile!(filename, read_step, species_list, reaction_list,
        system, speciesID)
    if (errcode == c_io_error) return errcode end
    print("End of input deck reading\n\n")

    return errcode
end


function ReadFile!(filename::String, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, speciesID::SpeciesID)
    
    errcode = 0

    open(filename,"r") do f
        line = 1
        while ! eof(f)
            s = readline(f)
            # ReadLine identifies each line on filename
            errcode = ReadLine!(s, read_step, species_list, reaction_list,
                system, speciesID)
            if (errcode == c_io_error)
                print("***ERROR*** Stop reading at file line ", line,"\n")
                return errcode  
            end
            line += 1
        end
    end
    return errcode
end


function ReadLine!(string::String, read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, speciesID::SpeciesID)

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
            errcode = StartBlock!(block_name, read_step, species_list,
                reaction_list, system, speciesID)
            if (errcode == c_io_error)
                print("***ERROR*** Something went wrong starting the ",
                    block_name, " block\n")
                return c_io_error
            end
        elseif (occursin("end", string))
            errcode = EndBlock!(block_name, read_step, species_list,
                reaction_list, system)
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
        errcode = ReadInputDeckEntry!(name, var, read_step,
            species_list, reaction_list, system, speciesID)
        if (errcode == c_io_error)
            print("***WARNING*** Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", name ," = ",var ,"\n")
        end
    elseif !(i_react === nothing)
        errcode = ReadInputDeckEntry!(string, string, read_step,
            species_list, reaction_list, system, speciesID) 
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


function StartFile!(read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, speciesID::SpeciesID)

    errcode = StartFile_Species!(read_step, species_list, speciesID) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input species block")
    end

    errcode = StartFile_Reactions!(reaction_list) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input reaction block")
    end
    
    errcode = StartFile_System!(read_step, system) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input system block")
    end
    
    return errcode
end


function StartBlock!(name::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, speciesID::SpeciesID)

    errcode = c_io_error
    if (occursin("system",name))
        global block_id = b_system
        errcode = StartSystemBlock!(read_step, system)
    elseif (occursin("species",name))
        global block_id = b_species
        errcode = StartSpeciesBlock!(read_step, species_list, speciesID)
    elseif (occursin("reactions",name))
        global block_id = b_reactions
        errcode = StartReactionsBlock!(reaction_list)
    end
    return errcode
end

function EndBlock!(name::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System)

    errcode = c_io_error
    global block_id = 0
    if (occursin("system",name))
        errcode = EndSystemBlock!(read_step, system)
    elseif (occursin("species",name))
        errcode = EndSpeciesBlock!(read_step, species_list)
    elseif (occursin("reactions",name))
        errcode = EndReactionsBlock!(reaction_list, species_list)
    end
    return errcode
end


function ReadInputDeckEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, speciesID::SpeciesID)

    errcode = c_io_error 
    
    if (block_id == b_system)
        errcode = ReadSystemEntry!(name, var, read_step, system)
    elseif (block_id == b_species)
        errcode = ReadSpeciesEntry!(name, var, read_step, species_list,
            speciesID)
    elseif (block_id == b_reactions)
        errcode = ReadReactionsEntry!(name, var, reaction_list,
            speciesID)
    end

    return errcode 
end

end