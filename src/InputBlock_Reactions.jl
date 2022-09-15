module InputBlock_Reactions 
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

using SharedData: Reaction, Species, SpeciesID, System
using SharedData: c_io_error
using SharedData: r_elastic, r_lower_threshold, r_diffusion, r_extended
using SharedData: r_emission_rate
using SharedData: K_to_eV, me, e 
using ParseReactions: LoadReaction!
using ParseReactions_PreRun: WriteRateCoefficientFunctions!
using PrintModule: PrintErrorMessage
using Printf: @printf, @sprintf


###############################################################################
################################  VARIABLES  ##################################
###############################################################################
global reaction_block_type = 0::Int64
global f_ReactionSet

function StartFile_Reactions!(read_step::Int64, reaction_list::Vector{Reaction}) 

    errcode = 0

    if read_step == 0
        global f_ReactionSet = open("src/ReactionSet.jl","w") 
        open("src/ReactionSet.Template","r") do f_temp
            while ! eof(f_temp)
                line_str = readline(f_temp, keep = true)
                if (line_str == "### START REACTION STRINGS ###\n")
                    break
                else
                    write(f_ReactionSet,line_str)
                end
            end
        end
    end

    return errcode
end


function StartReactionsBlock!(read_step::Int64, reaction_list::Vector{Reaction})

    errcode = 0
    global reaction_block_type = 0

    return errcode
end


function EndReactionsBlock!(read_step::Int64, reaction_list::Vector{Reaction},
    species_list::Vector{Species})

    errcode = 0 

    if read_step == 2
        # Set reacting neutral species
        for reaction in reaction_list
            errcode = IdentifyReactingNeutralSpecies!(reaction,
                species_list)
            if (errcode == c_io_error) return errcode end
        end
    end

    return errcode
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
        dens = Float64[]
        for s in species_list
            push!(temp, s.temp)
            push!(dens, s.dens)
        end

        for r in reaction_list
            ##### Test reaction charge balance
            charge_balance = 0.0
            i = 1
            for id in r.involved_species 
                fact = r.species_balance[i]
                charge_balance += species_list[id].charge * fact
                i += 1
            end
            if !(charge_balance == 0)
                error_str = @sprintf("Reaction %i: %s is not chaged balanced",r.id,r.name)
                PrintErrorMessage(system, error_str)
                return c_io_error
            end

            ##### Test elastic collisions
            if r.case == r_elastic
                if length(r.neutral_species_id) > 1
                    error_str = @sprintf("Elastic collision %i can only have one neutral reacting species", r.id)
                    PrintErrorMessage(system, error_str)
                    return c_io_error
                end
            end

            ##### Test reaction mass balance
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
                error_str = @sprintf("Reaction %i: %s is not mass balanced",r.id,r.name)
                PrintErrorMessage(system, error_str)
                return c_io_error
            end

            ##### Test rate coefficient
            for T_e in range(T_min, T_max, length=100)
                temp[sID.electron] = T_e 
                if system.prerun
                    if r.case == r_extended || r.case == r_diffusion || r.case == r_emission_rate
                        K = r.rate_coefficient(temp, dens, species_list, system, sID) 
                    else
                        K = r.rate_coefficient(temp, sID) 
                    end
                else
                    K = ReplaceExpressionValues(r.rate_coefficient, temp,
                        species_list, system, sID)
                end
                if K < 0.0
                    error_str = @sprintf("Reaction %i: %s with negative rate coeffciient at %.2f eV",r.id, r.name, T_e*K_to_eV)
                    PrintErrorMessage(system, error_str)
                    return c_io_error
                end
            end

            ##### Test energy threshold
            # - It is assumed that there is no E-threshold above +/-100 eV
            if abs(r.E_threshold) > 100.0/e
                error_str = @sprintf("Energy threshold is above 100 eV")
                PrintErrorMessage(system, error_str)
                return c_io_error
            end

            ##### Test emission reactions and self-absorption
            if r.case == r_emission_rate
                if length(r.reactant_species) > 1 || length(r.product_species) > 1
                    error_str = @sprintf("Radiation processes should only have one species at each side of the reaction")
                    PrintErrorMessage(system, error_str)
                    return c_io_error
                end

                if r.self_absorption
                    if r.g_high <= 0.0
                        error_str = @sprintf("Self absorption: statistical weight of higher state is not valid")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                    if r.g_low <= 0.0
                        error_str = @sprintf("Self absorption: statistical weight of lower state is not valid")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                    if r.g_high_total >= 1.e100
                        error_str = @sprintf("Self absorption: total statistical weight of higher states is not valid")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                    if r.g_low_total >= 1.e100
                        error_str = @sprintf("Self absorption: total statistical weight of lower states is not valid")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                    if r.wavelength <= 0.0
                        error_str = @sprintf("Self absorption: radiation wavelength is not valid")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                    if length(r.product_species) > 1
                        error_str = @sprintf("Self absorption: process only allows one product species")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                    if length(r.reactant_species) > 1
                        error_str = @sprintf("Self absorption: process only allows one reacting species")
                        PrintErrorMessage(system, error_str)
                        return c_io_error
                    end
                end
            end
        end

    elseif read_step == 0
        #  PRERUN
        write_flag = false 
        open("src/ReactionSet.Template","r") do f_temp
            while ! eof(f_temp)
                line_str = readline(f_temp, keep = true)

                if write_flag
                    write(f_ReactionSet,line_str)
                end

                if (line_str == "### END REACTION STRINGS ###\n")
                    write_flag = true
                end
                    
            end
        end
        close(f_ReactionSet)
    end

    return errcode
end


function ReadReactionsEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, reaction_list::Vector{Reaction}, system::System, 
    sID::SpeciesID)
    # Splits the input line contained in var into four parts
    # The fourth part is actually not necessary, but it is 
    # recommended to be included

    errcode = 0 

    if name == "reaction_type"
        if var == "elastic_scattering" || var == "elastic"
            global reaction_block_type = r_elastic
        elseif var == "diffusion"
            global reaction_block_type = r_diffusion
        elseif var == "lower_threshold"
            global reaction_block_type = r_lower_threshold
        elseif var == "emission_rate" || var == "emission"
            global reaction_block_type = r_emission_rate
        elseif var == "extended_function_parameters"
            global reaction_block_type = r_extended
        else
            errcode = c_io_error
            err_str = "Reaction block type not recognized" 
            if read_step > 1
                PrintErrorMessage(system, err_str)
            else
                print("***ERROR*** ", err_str,"\n")
            end
        end
        return errcode
    end

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

    # Second part: the threshold energy / emission parameters
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

    if read_step == 0 || read_step == 1

        # Initialize Reaction structure
        current_reaction = Reaction()
        InitializeReaction!(current_reaction, reaction_list)

        # Include description feature to the current reaction structure 
        if !(description_str === nothing)
            errcode = ParseDescription!(description_str, current_reaction)
            if (errcode == c_io_error) return errcode end
        end

        if read_step == 1
            # Main simulation procedure: set-up reaction structure
            errcode = LoadReaction!(current_reaction, system,
            reaction_process_str, e_threshold_str, sID)
            push!(reaction_list, current_reaction)
        elseif read_step == 0
            # PRERUN: read rate coefficients and write them to ReactionSet module
            errcode = WriteRateCoefficientFunctions!(current_reaction,
                rate_coeff_str, f_ReactionSet)
        end
    end

    return errcode
end


function InitializeReaction!(reaction::Reaction,
    reaction_list::Vector{Reaction})

    reaction.name = ""
    reaction.id = length(reaction_list) + 1 
    reaction.case = reaction_block_type
    reaction.neutral_species_id = Int64[]
    reaction.involved_species = Int64[]
    reaction.reactant_species = Int64[]
    reaction.product_species = Int64[] 
    reaction.species_balance= Int64[] 
    reaction.E_threshold = 0.0
    reaction.K_value = 0.0
    reaction.self_absorption = false
    reaction.g_high= 0.0
    reaction.g_low = 0.0
    reaction.g_high_total = 1.e100
    reaction.g_low_total = 1.e100
    reaction.wavelength = 0.0

end


function ParseDescription!(str::SubString{String}, reaction::Reaction)

    errcode = 0

    str = lowercase(str)

    if str == "elastic"
        reaction.case = r_elastic
    elseif str == "diffusion"
        reaction.case = r_diffusion
    elseif str == "lower_threshold"
        reaction.case = r_lower_threshold
    elseif str == "extended_function_parameters"
        reaction.case = r_extended
    elseif str == ""
        errcode = 0
    else
        errcode = c_io_error
        print("***ERROR*** Reaction description '",str,"' was not recognized\n")
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

end