module InputData

using SharedData: c_io_error
using SharedData: p_icp_id, p_ccp_id 
using SharedData: Species, Reaction, System, SpeciesID

using InputBlock_Species: StartFile_Species!
using InputBlock_Species: StartSpeciesBlock!, EndSpeciesBlock!, ReadSpeciesEntry!

using InputBlock_Reactions: StartFile_Reactions!
using InputBlock_Reactions: StartReactionsBlock!, EndReactionsBlock!, ReadReactionsEntry!
using InputBlock_Reactions: r_elastic, r_wall_loss 

using InputBlock_System: StartFile_System!
using InputBlock_System: StartSystemBlock!, EndSystemBlock!, ReadSystemEntry!

using PlasmaParameters: GetLambda, GetGamma
using EvaluateExpressions: ReplaceSystemSymbolS!

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

function SetupInputData!(filename::String, setup_flag::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, speciesID::SpeciesID)

    errcode = c_io_error

    errcode = ReadInputData!(filename, setup_flag,
        species_list, reaction_list, system, speciesID)
    if (errcode == c_io_error)
        print("***ERROR*** Failed to read the input deck. Abort code\n")
    end

    return errcode
end


function ReadInputData!(filename::String, setup_flag::Int64,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, speciesID::SpeciesID)

    errcode = c_io_error 

    if (setup_flag == setup_pre_run)
        read_step = 0 
        print("Reading $read_step of the input deck...\n")
        errcode = StartFile!(read_step, species_list, reaction_list, system,
            speciesID)
        if (errcode == c_io_error) return errcode end
        errcode = ReadFile!(filename, read_step, species_list, reaction_list,
            system, speciesID)
        if (errcode == c_io_error) return errcode end
        print("End of input deck reading\n\n")

    elseif (setup_flag == setup_main_run)
        # Opens the file given in filename and reads each line
        for read_step in 1:2
            print("Reading $read_step of the input deck...\n")
            errcode = StartFile!(read_step, species_list, reaction_list,
                system, speciesID)
            if (errcode == c_io_error) return errcode end
            errcode = ReadFile!(filename, read_step, species_list,
                reaction_list, system, speciesID)
            if (errcode == c_io_error) return errcode end
            print("End of input deck reading\n\n")
        end
        errcode = CheckSpeciesList(species_list, reaction_list, system)
        errcode = CheckReactionList(species_list, reaction_list, system)
    else
        print("***ERROR*** Setup flag is not recognized\n")
    end
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

    errcode = StartFile_Reactions!(read_step, reaction_list) 
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
        errcode = StartReactionsBlock!(read_step, reaction_list)
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
        errcode = EndReactionsBlock!(read_step, reaction_list, species_list)
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
        errcode = ReadReactionsEntry!(name, var, read_step, reaction_list,
            speciesID)
    end

    return errcode 
end


function CheckSpeciesList(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

    errcode = 0

    for s in species_list
        s_id = s.id
        for r in reaction_list
            # is species s involved?
            i_involved = findall(x->x==s_id, r.reactant_species)
            if i_involved==Int64[]
                continue
            else
                push!(s.reaction_list, r)
            end
        end
        if system.power_input_method == p_icp_id
            s.Lambda = GetLambda(system)
            s.gamma = GetGamma()
        end
    end
    return errcode
end


function CheckReactionList(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

    errcode = 0

    for r in reaction_list
        # Update rate coefficient values
        if !(typeof(r.rate_coefficient) == Float64)
            ReplaceSystemSymbolS!(r.rate_coefficient, system)
        end

        # Test reaction charge balance
        if !(r.case == r_wall_loss)
            charge_balance = 0.0
            i = 1
            for id in r.involved_species 
                fact = r.species_balance[i]
                charge_balance += species_list[id].charge * fact
                i += 1
            end
            if !(charge_balance == 0)
                print("***ERROR*** Reaction ",r.id," is unbalanced\n")
                return c_io_error
            end
        end

        # Test elastic collisions
        if r.case == r_elastic
            if length(r.neutral_species_id) > 1
                print("***ERROR*** Elastic collision ",r.id,
                    " can only have one neutral reacting species\n")
                return c_io_error
            end
        end
    end
    return errcode
end

end