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

using SharedData: c_io_error, e
using SharedData: Species, Reaction, System, SpeciesID
using SharedData: r_diffusion, r_emission_rate
using ReactionSet: K_funct_list
using EvaluateExpressions: ReplaceConstantValues!
using PrintModule: PrintErrorMessage 



function LoadReaction!(current_reaction::Reaction, system::System,
    reaction_process_str::SubString{String},
    e_threshold_str::SubString{String}, sID::SpeciesID)

    errcode = 0

    # Parse each term of the input line
    current_reaction.name = reaction_process_str

    errcode = ParseReaction!(reaction_process_str, current_reaction, sID)
    if (errcode == c_io_error) return errcode end

    errcode = ParseEThreshold!(e_threshold_str, current_reaction)
    if (errcode == c_io_error) return errcode end

    if system.prerun
        #if current_reaction.case == r_diffusion
        #    current_reaction.rate_coefficient =
        #elseif current_reaction.case == r_emission_rate
        #    current_reaction.rate_coefficient = 
        #else
            current_reaction.rate_coefficient = K_funct_list[current_reaction.id]
        #end
    else
        errcode = GetRateCoefficientExpr(rate_coeff_str, current_reaction)
        if (errcode == c_io_error) return errcode end
    end

    return errcode
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

end