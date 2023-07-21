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

module ParseReactions

using SharedData: c_io_error, r_emission_rate
using SharedData: e
using SharedData: Species, Reaction, System, SpeciesID
using ReactionSet: K_funct_list
using EvaluateExpressions: ReplaceConstantValues!
using Printf: @sprintf
using PrintModule: PrintErrorMessage 
using InputBlock_System: GetUnits!



function LoadReaction!(current_reaction::Reaction, system::System,
    reaction_process_str::SubString{String},
    secondentry_str::SubString{String}, sID::SpeciesID)

    errcode = 0

    # Parse each term of the input line
    current_reaction.name = reaction_process_str

    errcode = ParseReaction!(reaction_process_str, current_reaction, system,
        sID)
    if (errcode == c_io_error) return errcode end

    errcode = ParseSecondEntry!(secondentry_str, current_reaction, system)
    if (errcode == c_io_error) return errcode end

    if system.prerun
        current_reaction.rate_coefficient = K_funct_list[current_reaction.id]
    else
        errcode = GetRateCoefficientExpr(rate_coeff_str, current_reaction, system)
        if (errcode == c_io_error) return errcode end
    end

    return errcode
end


function ParseReaction!(str::SubString{String}, reaction::Reaction,
    system::System, speciesID::SpeciesID)

    errcode = c_io_error

    idx = findfirst("->",str)
    if !(idx===nothing)
        # Get reactant and product strings
        reactant_str = strip(str[1:idx[1]-1])
        product_str = strip(str[idx[2]+1:end])
        
        # Get Vector{Int64} with species IDs
        input_rea_s = GetSpeciesFromString!(reactant_str, speciesID)
        if (input_rea_s == c_io_error)
            err_message = "While getting reacting species"
            PrintErrorMessage(system, err_message)
            return errcode
        end

        input_pro_s = GetSpeciesFromString!(product_str, speciesID)
        if (input_pro_s == c_io_error)
            err_message = "While getting product species"
            PrintErrorMessage(system, err_message)
            return errcode
        end

        # Get balance, reactant and involved species vectors
        errcode = GetReactionSpeciesLists!(input_rea_s, input_pro_s, reaction,
            system)
        if (input_pro_s == c_io_error)
            err_message = "While getting species reaction balance, reactants and species vector"
            PrintErrorMessage(system, err_message)
            return errcode
        end
    end
    return errcode
end


function ParseSecondEntry!(str::SubString{String}, reaction::Reaction,
    system::System)

    # Must provide: E_threshold
    errcode = 0

    try
        if reaction.case == r_emission_rate
            # String should be of the form "[X,Y,Z]" where X,Y,Z are numbers

            # Check whether there is a vector
            left_braket = findfirst("[",str)
            if left_braket === nothing
                units, new_str = GetUnits!(str)
                reaction.E_threshold = parse(Float64, new_str) * units
            else
                reaction.self_absorption = true
                GetEmissionParameters!(reaction, str)
            end
        else
            units, new_str = GetUnits!(str)
            reaction.E_threshold = parse(Float64, new_str) * units
        end
    catch
        if reaction.case == r_emission_rate
            err_message = "While parsing statistical weights and wavelength"
        else
            err_message = "While parsing E threshold"
        end
        PrintErrorMessage(system, err_message)
        errcode = c_io_error
    end
        
    return errcode
end


function GetEmissionParameters!(reaction::Reaction, str::SubString{String})

    # Identify vector brakets
    left_braket = findfirst("[",str)
    left_ix = left_braket[1] + 1
    right_braket = findlast("]",str)
    right_ix = right_braket[1] - 1
    str_no_brakets = str[left_ix:right_ix]

    # Find commas
    comma_ix_list = Int64[]
    start_ix = 1
    for i in range(1,4,step=1)
        comma = findnext(",", str_no_brakets,start_ix)
        ix = comma[1]
        push!(comma_ix_list, ix)
        start_ix = ix+1
    end
    c_1 = comma_ix_list[1]
    c_2 = comma_ix_list[2]
    c_3 = comma_ix_list[3]
    c_4 = comma_ix_list[4]

    # Identify parameters
    g_p_str = str_no_brakets[1:c_1-1]
    units, g_p_str = GetUnits!(g_p_str)
    reaction.g_high = parse(Float64, g_p_str) * units

    g_k_str = str_no_brakets[c_1+1:c_2-1]
    units, g_k_str = GetUnits!(g_k_str)
    reaction.g_low = parse(Float64, g_k_str) * units

    g_p_tot_str = str_no_brakets[c_2+1:c_3-1]
    units, g_p_tot_str = GetUnits!(g_p_tot_str)
    reaction.g_high_total = parse(Float64, g_p_tot_str) * units

    g_k_tot_str = str_no_brakets[c_3+1:c_4-1]
    units, g_k_tot_str = GetUnits!(g_k_tot_str)
    reaction.g_low_total = parse(Float64, g_k_tot_str) * units

    wavelen_str = str_no_brakets[c_4+1:end]
    units, wavelen_str = GetUnits!(wavelen_str)
    reaction.wavelength = parse(Float64,wavelen_str) * units

end


function GetSpeciesFactor!(str::SubString{String})

    # This functions checks wether each species in the
    # reaction has a multiple on its left

    fact = 1
    num_list = range(2,9,step=1)
    for num in num_list 
        num_str = string(num)
        if str[1:1] == num_str 
            fact = num
            str = str[2:end]
            break
        end
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
    reaction::Reaction, system::System)

    errcode = 0
    try
        # Loop over species list and
        #  - if not in involved_species -> push
        #  - if already in, move on
        # Loop for reactant/product species
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
            end
            push!(reaction.reactant_species, r)
        end
        
        # Loop for product species
        for p in prod 
            # Check involved species
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

            # Check product species
            already_in = false
            i = 0
            for is in reaction.product_species
                i += 1
                if p == is
                    already_in = true
                    break
                end
            end
            if !already_in
                push!(reaction.product_species, p)
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
        PrintErrorMessage(system, "While setting up reaction species lists")
        errcode = c_io_error
    end
    return errcode
end


function GetRateCoefficientExpr(str::SubString{String},
    reaction::Reaction, system::System)

    errcode = 0 

    try
        expr = Meta.parse(str)
        if typeof(expr)==Expr
            ReplaceConstantValues!(expr)
        end

        reaction.rate_coefficient = expr

    catch
        errcode = c_io_error
        PrintErrorMessage(system, "While getting rate coefficient expression")
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
                err_str = @sprintf("Reaction species %s is not recognized", str)
                PrintErrorMessage(system, err_str)
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
                err_str = @sprintf("Reaction species %s is not recognized", s)
                PrintErrorMessage(system, err_str)
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

end