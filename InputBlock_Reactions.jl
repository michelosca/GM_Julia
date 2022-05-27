# Copyright (C) 2021 Michel Osca Engelbrecht
#
# This file is part of GM Julia.
#
# GM Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GM Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module InputBlock_Reactions

using SharedData: c_io_error, e, me, K_to_eV
using SharedData: Species, Reaction, System, SpeciesID
using SharedData: r_elastic, r_wall_loss
using ReactionSet: K_funct_list
using EvaluateExpressions: ReplaceConstantValues!, ReplaceSystemSymbols!
using EvaluateExpressions: ReplaceSpeciesSymbols!, ReplaceTempSymbols!
using EvaluateExpressions: ReplaceDensSymbols!
using EvaluateExpressions: ReplaceExpressionValues
using Printf



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
    read_step::Int64, reaction_list::Vector{Reaction}, system::System, 
    speciesID::SpeciesID)
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
        
        if system.prerun
            current_reaction.rate_coefficient = K_funct_list[current_reaction.id]
        else
            errcode = GetRateCoefficientExpr(rate_coeff_str, current_reaction)
            if (errcode == c_io_error) return errcode end
        end

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
    elseif str[1:1] == "3"
        fact = 3
        str = str[2:end]
    elseif str[1:1] == "4"
        fact = 4
        str = str[2:end]
    elseif str[1:1] == "5"
        fact = 5
        str = str[2:end]
    end
    return fact, str
end


function SelectSpeciesID(s::SubString{String}, speciesID::SpeciesID)

    id = 0
    # ARGON 
    if (s == "Ar")
        id = speciesID.Ar 
    elseif (s == "Ar+")
        id = speciesID.Ar_Ion 
    elseif (s == "Ar_m")
        id = speciesID.Ar_m 
    elseif (s == "Ar_r")
        id = speciesID.Ar_r 
    elseif (s == "Ar_4p")
        id = speciesID.Ar_4p 

    # ATOMIC OXYGEN 
    elseif (s == "O")
        id = speciesID.O 
    elseif (s == "O+")
        id = speciesID.O_Ion 
    elseif (s == "O-")
        id = speciesID.O_negIon 
    elseif (s == "O_1s")
        id = speciesID.O_1s
    elseif (s == "O_3s")
        id = speciesID.O_3s
    elseif (s == "O_5s")
        id = speciesID.O_5s
    elseif (s == "O_1d")
        id = speciesID.O_1d
    elseif (s == "O_3p")
        id = speciesID.O_3p
    elseif (s == "O_5p")
        id = speciesID.O_5p

    # MOLECULAR OXYGEN
    elseif (s == "O2")
        id = speciesID.O2
    elseif (s == "O2_v")
        id = speciesID.O2_v
    elseif (s == "O2+")
        id = speciesID.O2_Ion 
    elseif (s == "O2-")
        id = speciesID.O2_negIon 
    elseif (s == "O2_a1Ag")
        id = speciesID.O2_a1Ag 
    elseif (s == "O2_a1Ag_v")
        id = speciesID.O2_a1Ag_v 
    elseif (s == "O2_b1Su")
        id = speciesID.O2_b1Su 
    elseif (s == "O2_b1Su_v")
        id = speciesID.O2_b1Su_v 

    # OZONE and O4 
    elseif (s == "O3")
        id = speciesID.O3
    elseif (s == "O3_v")
        id = speciesID.O3_v
    elseif (s == "O3+")
        id = speciesID.O3_Ion
    elseif (s == "O3-")
        id = speciesID.O3_negIon
    elseif (s == "O4+")
        id = speciesID.O4_Ion
    elseif (s == "O4-")
        id = speciesID.O4_negIon

    # ELECTRON
    elseif (s == "e")
        id = speciesID.electron 
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


function GetRateCoefficientExpr(str::SubString{String},
    reaction::Reaction)

    errcode = 0 

    try
        expr = Meta.parse(str)
        if typeof(expr)==Expr
            ReplaceConstantValues!(expr)
        end

        reaction.rate_coefficient = expr

    catch
        errcode = c_io_error
        print("***ERROR*** While getting rate coefficient expression\n")
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
    species_list::Vector{Species}, system::System, sID::SpeciesID)

    errcode = 0
    if read_step == 2
        print("Test reactions...\n")
        print("  - Charge balance\n")
        print("  - Mass balance\n")
        print("  - Elastic scattering reactions\n")
        T_min = 1.5 / K_to_eV
        T_max = 4.0 / K_to_eV
        @printf("  - Rate coefficients values between %.2f <= T_e <= %.2f [eV]\n",
            T_min*K_to_eV, T_max*K_to_eV)
        # Setup temp array, used later when testing rate coefficients
        temp = Float64[]
        for s in species_list
            push!(temp, s.temp)
        end

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
                    print("***ERROR*** Reaction ",r.id,": ",r.name," is charged unbalanced\n")
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

            # Test reaction mass balance
            mass_balance = 0.0
            i = 1
            for id in r.involved_species 
                if id == sID.electron
                    i += 1
                    continue
                end
                fact = r.species_balance[i]
                mass_balance += species_list[id].mass * fact
                i += 1
            end
            if !(abs(mass_balance) <= me)
                print("***ERROR*** Reaction ",r.id,": ",r.name," is mass unbalanced\n")
                return c_io_error
            end

            # Test rate coefficient
            for T_e in range(T_min, T_max, length=100)
                temp[sID.electron] = T_e 
                if system.prerun
                    if r.case == r_wall_loss
                        K = r.rate_coefficient(temp, species_list, system, sID) 
                    else
                        K = r.rate_coefficient(temp, sID) 
                    end
                else
                    K = ReplaceExpressionValues(r.rate_coefficient, temp,
                        species_list, system, sID)
                end
                if K < 0.0
                    @printf("***ERROR*** Reaction %i: %s has negative rate coefficient value at %.2f eV\n",r.id, r.name, T_e * K_to_eV)
                    return c_io_error
                end
            end
        end
    end

    return errcode
end

end