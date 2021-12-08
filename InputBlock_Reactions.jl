module InputBlock_Reactions

using SharedData: c_io_error, e
using SharedData: Species, Reaction, SpeciesID
using SharedData: r_energy_sink, r_elastic, r_wall_loss
using ReactionSet: K_funct_list


###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - StartFile_Reactions
# - StartReactionsBlock
# - ReadReactionsEntry
#   - ParseReaction
#     - GetSpeciesFromString
#       - SelectSpeciesID
#     - GetReactionSpeciesLists
#   - ParseEThreshold
#   - ParseRateCoefficient
#   - ParseDescription
#   - AddReactionToList 
# - EndReactionsBlock

function StartFile_Reactions!(read_step::Int64, reaction_list::Vector{Reaction}) 

    errcode = 0

    return errcode
end

function StartReactionsBlock!(read_step::Int64, reaction_list::Vector{Reaction})

    errcode = 0

    return errcode
end


function ReadReactionsEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, reaction_list::Vector{Reaction}, speciesID::SpeciesID)
    # Splits the input line contained in var into four parts
    # The fourth part is actually not necessary, but it is 
    # recommended to be included

    errcode = 0 

    # First part: the reaction process
    idx = findfirst(";", var)
    if !(idx === nothing)
        idx = idx[1]
        reaction_process_str = strip(var[1:idx-1])
        var = var[idx+1:end]
        str_track = true
    else
        errcode = c_io_error
        str_track = false 
    end

    # Second part: the threshold energy
    idx = findfirst(";", var)
    if (!(idx === nothing) && str_track)
        idx = idx[1]
        e_threshold_str = strip(var[1:idx-1])
        var = var[idx+1:end]
        str_track = true
    else
        errcode = c_io_error
        str_track = false 
    end

    # Third part: the rate coefficient expression
    idx = findfirst(";", var)
    if (!(idx === nothing) && str_track)
        idx = idx[1]
        rate_coeff_str = strip(var[1:idx-1])
        var = var[idx+1:end]
        str_track = true
    else
        errcode = c_io_error
        str_track = false
    end

    # Fourth part, if existing: reaction description string
    if (str_track)
        description_str = strip(var)
    else
        description_str = nothing 
    end
    if (errcode == c_io_error) return errcode end

    if (read_step ==1)

        # Parse each term of the input line
        current_reaction = Reaction()
        InitializeReaction!(current_reaction, reaction_list)
        current_reaction.name = reaction_process_str

        errcode = ParseReaction!(reaction_process_str, current_reaction,
            speciesID)
        if (errcode == c_io_error) return errcode end

        errcode = ParseEThreshold!(e_threshold_str, current_reaction)
        if (errcode == c_io_error) return errcode end

        if !(description_str === nothing)
            errcode = ParseDescription!(description_str, current_reaction)
            if (errcode == c_io_error) return errcode end
        end
        
        current_reaction.rate_coefficient = K_funct_list[current_reaction.id]

        # Add current_reaction to reaction_list
        push!(reaction_list, current_reaction)
    end

    return errcode
end


function InitializeReaction!(reaction::Reaction, reaction_list::Vector{Reaction})

    reaction.name = ""
    reaction.id = length(reaction_list) + 1 
    reaction.case = 0
    reaction.neutral_species_id = Int64[]
    reaction.E_threshold = 0.0

end


function ParseReaction!(str::SubString{String}, reaction::Reaction,
    speciesID::SpeciesID)

    errcode = c_io_error

    idx = findfirst("->",str)
    if !(idx===nothing)
        # Get reactant and product strings
        reactant_str = strip(str[1:idx[1]-1])
        product_str = strip(str[idx[2]+1:end])
        
        # Get Vector{Int64} with species IDs
        input_rea_s = GetSpeciesFromString!(reactant_str, speciesID)
        if (input_rea_s == c_io_error) return errcode end

        input_pro_s = GetSpeciesFromString!(product_str, speciesID)
        if (input_pro_s == c_io_error) return errcode end

        # Get balance, reactant and involved species vectors
        errcode = GetReactionSpeciesLists!(input_rea_s, input_pro_s, reaction)
        if (input_pro_s == c_io_error) return errcode end
    end
    return errcode
end


function GetSpeciesFromString!(str::SubString{String}, speciesID::SpeciesID)
    # This function links the species found in the given string (str)
    # to the species ids predefined in the code

    s_list = Int64[]
    next_species = true 
    while next_species
        idx = findfirst(" + ", str)
        if (idx===nothing)
            fact, str= GetSpeciesFactor!(str)
            s_id = SelectSpeciesID(str, speciesID)
            if (s_id == 0)
                print("***ERROR*** Reaction species ",str ," is not recognized\n")
                return c_io_error
            else 
                while fact >= 1
                    push!(s_list, s_id)
                    fact -= 1
                end
            end
            next_species = false
        else
            s = strip(str[1:idx[1]-1])
            str = strip(str[idx[2]+1:end])
            fact, s = GetSpeciesFactor!(s)
            s_id = SelectSpeciesID(s, speciesID)
            if (s_id == 0)
                print("***ERROR*** Reaction species ",s ," is not recognized\n")
                return c_io_error
            else 
                while fact >= 1
                    push!(s_list, s_id)
                    fact -= 1
                end
            end
        end
    end
    return s_list
end


function GetSpeciesFactor!(str::SubString{String})

    fact = 1
    if str[1:1] == "2"
        fact = 2
        str = str[2:end]
    end
    return fact, str
end


function SelectSpeciesID(s::SubString{String}, speciesID::SpeciesID)

    id = 0
    # Neutral species
    if (s == "Ar")
        id = speciesID.Ar 
    elseif (s == "O")
        id = speciesID.O 
    elseif (s == "O2")
        id = speciesID.O2

    # Charged and excited species
    elseif (s == "Ar+")
        id = speciesID.Ar_Ion 
    elseif (s == "Ar_m")
        id = speciesID.Ar_m 
    elseif (s == "Ar_r")
        id = speciesID.Ar_r 
    elseif (s == "Ar_4p")
        id = speciesID.Ar_4p 
    elseif (s == "e")
        id = speciesID.electron 
    elseif (s == "O+")
        id = speciesID.O_Ion 
    elseif (s == "O-")
        id = speciesID.O_negIon 
    elseif (s == "O2+")
        id = speciesID.O2_Ion 
    elseif (s == "O_1d")
        id = speciesID.O_1d 
    elseif (s == "O2_a1Ag")
        id = speciesID.O2_a1Ag 
    end
    return id 
end


function GetReactionSpeciesLists!(reac::Vector{Int64}, prod::Vector{Int64},
    reaction::Reaction)

    errcode = 0
    try
        reaction.involved_species = Int64[]
        reaction.reactant_species = Int64[]

        # Loop over species list and
        #  - if not in involved_species -> push
        #  - if already in, move on
        # Loop for reactant species
        for r in reac
            already_in = false
            i = 0
            for is in reaction.involved_species
                i += 1
                if r == is
                    already_in = true
                    break
                end
            end
            if !already_in
                push!(reaction.involved_species, r)
                push!(reaction.reactant_species, r)
            end
        end
        
        # Loop for product species
        for p in prod 
            already_in = false
            i = 0
            for is in reaction.involved_species
                i += 1
                if p == is
                    already_in = true
                    break
                end
            end
            if !already_in
                push!(reaction.involved_species, p)
            end
        end

        # Get species balance between reactants and products
        n_species = length(reaction.involved_species)
        reaction.species_balance= zeros(Int64, n_species) 
        for i in 1:n_species
            s = reaction.involved_species[i]
            for r in reac
                if s==r
                    reaction.species_balance[i] -= 1
                end
            end

            for p in prod
                if s==p
                    reaction.species_balance[i] += 1
                end
            end
        end

    catch
        print("***ERROR*** While setting up reaction species lists\n")
        errcode = c_io_error
    end
    return errcode
end


function ParseEThreshold!(str::SubString{String}, reaction::Reaction)

    # Must provide: E_threshold
    errcode = 0

    try
        # Default energy units: eV
        units_fact = e
        if (occursin("J", str))
            units_fact = 1.0
            idx = findfirst("J", str)
            idx = idx[1]-1
            str = string(str[1:idx])
        elseif (occursin("eV", str))
            idx = findfirst("eV", str)
            idx = idx[1]-1
            str = string(str[1:idx])
        end
        reaction.E_threshold = parse(Float64, str) * units_fact
    catch
        print("***ERROR*** While parsing E threshold\n")
        errcode = c_io_error
    end
        
    return errcode
end


function ParseDescription!(str::SubString{String}, reaction::Reaction)

    errcode = 0
    reaction.case = 0

    str = lowercase(str)

    if (str == "elastic")
        reaction.case = r_elastic
    elseif (str == "energy_sink" || str == "energy sink")
        reaction.case = r_energy_sink
    elseif (str == "wall_rate_coefficient")
        reaction.case = r_wall_loss
    elseif (str == "")
        errcode = 0
    else
        errcode = c_io_error
        print("***ERROR*** Reaction description was not recognized\n")
    end

    return errcode
end


function EndReactionsBlock!(read_step::Int64, reaction_list::Vector{Reaction},
    species_list::Vector{Species})
    errcode = 0 

    if (read_step == 2)
        # Set reacting neutral species
        for reaction in reaction_list
            errcode = IdentifyReactingNeutralSpecies!(reaction,
                species_list)
            if (errcode == c_io_error) return errcode end
        end
    end

    return errcode
end


function IdentifyReactingNeutralSpecies!(reaction::Reaction,
    species_list::Vector{Species})

    for react_id in reaction.reactant_species
        species = species_list[react_id]
        if species.charge == 0
            already_in = false
            for n_id in reaction.neutral_species_id
                if react_id == n_id
                    already_in = true
                end
            end
            if !already_in
                push!(reaction.neutral_species_id, react_id)
            end
        end
    end
end


function EndFile_Reactions!(read_step::Int64, reaction_list::Vector{Reaction},
    species_list::Vector{Species})

    errcode = 0
    if read_step == 2
        for r in reaction_list
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
    end

    return errcode
end

end