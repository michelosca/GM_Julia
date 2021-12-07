module InputData

using SharedData: c_io_error
using SharedData: Species, Reaction, System, SpeciesID, OutputBlock

using InputBlock_Species: StartFile_Species!, EndFile_Species!
using InputBlock_Species: StartSpeciesBlock!, EndSpeciesBlock!, ReadSpeciesEntry!

using InputBlock_Reactions: StartFile_Reactions!, EndFile_Reactions!
using InputBlock_Reactions: StartReactionsBlock!, EndReactionsBlock!, ReadReactionsEntry!
using InputBlock_Reactions: r_elastic, r_wall_loss 

using InputBlock_System: StartFile_System!, EndFile_System!
using InputBlock_System: StartSystemBlock!, EndSystemBlock!, ReadSystemEntry!

using InputBlock_Output: StartFile_Output!, EndFile_Output!
using InputBlock_Output: StartOutputBlock!, EndOutputBlock!, ReadOutputEntry!

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
# INPUT BLOCK IDs
global block_id = 0
const b_system = 1
const b_species = 2
const b_reactions = 3
const b_output = 4
const b_constants = 5

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

function SetupInputData!(filename::String,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, speciesID::SpeciesID)

    errcode = c_io_error

    errcode = ReadInputData!(filename,
        species_list, reaction_list, system, output_list, speciesID)
    if (errcode == c_io_error)
        print("***ERROR*** Failed to read the input deck. Abort code\n")
    end

    return errcode
end


function ReadInputData!(filename::String,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, speciesID::SpeciesID)

    errcode = 0

    # Opens the file given in filename and reads each line
    for read_step in 1:2
        print("Reading $read_step of the input deck...\n")
        errcode = StartFile!(read_step, species_list, reaction_list,
            system, output_list, speciesID)
        if (errcode == c_io_error) return errcode end

        errcode = ReadFile!(filename, read_step, species_list,
            reaction_list, system, output_list, speciesID)
        if (errcode == c_io_error) return errcode end

        errcode = EndFile!(read_step, species_list, reaction_list,
            system, output_list, speciesID)
        if (errcode == c_io_error) return errcode end
        print("End of input deck reading\n\n")
    end

    return errcode
end


function ReadFile!(filename::String, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, speciesID::SpeciesID)
    
    errcode = 0

    open(filename,"r") do f
        constants = Tuple{SubString{String}, SubString{String}}[]
        line = 1
        while ! eof(f)
            s = readline(f)
            # ReadLine identifies each line on filename
            errcode = ReadLine!(s, read_step, species_list, reaction_list,
                system, output_list, speciesID, constants)
            if (errcode == c_io_error)
                print("***ERROR*** Stop reading at file line ", line,"\n")
                return errcode  
            end
            line += 1
        end
    end
    return errcode
end


function ReadLine!(str::String, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, speciesID::SpeciesID,
    constants::Vector{Tuple{SubString{String},SubString{String}}})

    errcode = c_io_error 

    # Trimm out commets
    i_comment = findfirst("#", str)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        str = str[begin:i_comment-1]
    end
    str = strip(str)

    # Signal wether this is a begin/end:block 
    i_block = findfirst(":", str)
    # Signla whether it is a "name = var" line
    i_eq = findfirst("=", str)
    # Signla whether it is a reaction line
    i_react = findfirst(";", str)

    # Check line
    if !(i_block === nothing)
        i_block = i_block[1]
        global block_name = str[i_block+1:end]
        if (occursin("begin", str))
            errcode = StartBlock!(block_name, read_step, species_list,
                reaction_list, system, output_list, speciesID)
            if (errcode == c_io_error)
                print("***ERROR*** Something went wrong starting the ",
                    block_name, " block\n")
                return c_io_error
            end
        elseif (occursin("end", str))
            errcode = EndBlock!(block_name, read_step, species_list,
                reaction_list, system, output_list, speciesID)
            if (errcode == c_io_error)
                print("***ERROR*** Something went wrong ending the ",
                    block_name, " block\n")
                return c_io_error
            end
        end
    elseif !(i_eq === nothing)
        i_eq = i_eq[1]
        name = strip(str[begin:i_eq-1])
        var = strip(str[i_eq+1:end])
        errcode = ReadInputDeckEntry!(name, var, read_step,
            species_list, reaction_list, system, output_list, speciesID,
            constants)
        if (errcode == c_io_error)
            print("***WARNING*** Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", name ," = ",var ,"\n")
        end
    elseif !(i_react === nothing)
        errcode = ReadInputDeckEntry!(str, str, read_step,
            species_list, reaction_list, system, output_list, speciesID,
            constants) 
        if (errcode == c_io_error)
            print("***WARNING*** Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", str ,"\n")
        end
    else
        errcode = 0
    end
    return errcode 
end


function StartFile!(read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System,
    output_list::Vector{OutputBlock}, speciesID::SpeciesID)

    errcode = StartFile_Species!(read_step, species_list, speciesID) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input species block")
    end

    errcode = StartFile_Reactions!(read_step, reaction_list) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input reaction block")
    end
    
    errcode = StartFile_System!(read_step, system) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input system block")
    end
    
    errcode = StartFile_Output!(read_step, output_list) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input output block")
    end
    
    return errcode
end


function EndFile!(read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System,
    output_list::Vector{OutputBlock}, speciesID::SpeciesID)

    errcode = EndFile_Species!(read_step, species_list, reaction_list, system,
        speciesID)
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input species block")
    end

    errcode = EndFile_Reactions!(read_step, reaction_list, species_list)
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input reaction block")
    end
    
    errcode = EndFile_System!(read_step, system) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input system block")
    end
    
    errcode = EndFile_Output!(read_step, output_list, species_list,
        reaction_list) 
    if (errcode == c_io_error)
        print("***ERROR*** While initializing the input output block")
    end
    
    return errcode
end


function StartBlock!(name::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, speciesID::SpeciesID)

    errcode = c_io_error
    if (occursin("system",name))
        global block_id = b_system
        errcode = StartSystemBlock!(read_step, system)
    elseif (occursin("species",name))
        global block_id = b_species
        errcode = StartSpeciesBlock!(read_step, species_list, speciesID)
    elseif (occursin("reactions",name))
        global block_id = b_reactions
        errcode = StartReactionsBlock!(read_step, reaction_list)
    elseif (occursin("output",name))
        global block_id = b_output
        errcode = StartOutputBlock!(read_step, output_list)
    elseif (occursin("constants",name))
        global block_id = b_constants
        errcode = 0
    end
    return errcode
end

function EndBlock!(name::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, sID::SpeciesID)

    errcode = c_io_error
    global block_id = 0
    if (occursin("system",name))
        errcode = EndSystemBlock!(read_step, system)
    elseif (occursin("species",name))
        errcode = EndSpeciesBlock!(read_step, species_list, sID)
    elseif (occursin("reactions",name))
        errcode = EndReactionsBlock!(read_step, reaction_list, species_list)
    elseif (occursin("output",name))
        errcode = EndOutputBlock!(read_step, output_list)
    elseif (occursin("constants",name))
        errcode = 0
    end
    return errcode
end


function ReadInputDeckEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System,
    output_list::Vector{OutputBlock}, speciesID::SpeciesID,
    constants::Vector{Tuple{SubString{String},SubString{String}}})

    errcode = c_io_error 

    var = CheckConstantValues!(var, constants)
    
    if (block_id == b_system)
        errcode = ReadSystemEntry!(name, var, read_step, system)
    elseif (block_id == b_species)
        errcode = ReadSpeciesEntry!(name, var, read_step, species_list,
            speciesID)
    elseif (block_id == b_reactions)
        errcode = ReadReactionsEntry!(name, var, read_step, reaction_list,
            speciesID)
    elseif (block_id == b_output)
        errcode = ReadOutputEntry!(name, var, read_step, output_list,
            species_list)
    elseif (block_id == b_constants)
        push!(constants, (name,var))
        errcode = 0
    end

    return errcode 
end


function CheckConstantValues!(var::SubString{String}, constants::Vector{Tuple{SubString{String},SubString{String}}})

    for c in constants
        if var == c[1]
            var = c[2]
            break
        end
    end

    return var
end

end