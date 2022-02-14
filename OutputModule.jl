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

module OutputModule

using SharedData: Species, Reaction, System, SpeciesID, OutputBlock
using SharedData: kb, K_to_eV
using SharedData: c_io_error, r_wall_loss
using SharedData: o_scale_lin, o_scale_log
using SharedData: o_single_run, o_pL, o_dens, o_temp, o_power, o_pressure
using SharedData: o_pressure_percent, neutral_species_id
using EvaluateExpressions: ReplaceExpressionValues
using SolveSystem: ExecuteProblem
using Printf
using PrintModule: PrintSpeciesList
###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GenerateOutputs
#   - SetupInitialConditions
#   - Output_PL

function GenerateOutputs!(
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, sID::SpeciesID)
    
    errcode = 0

    for output in output_list

        if output.case[1] == o_single_run
            sol = @time ExecuteProblem(species_list, reaction_list, system, sID)
            errcode = @time LoadOutputBlock!(output, sol, species_list,
                reaction_list, system, sID)
            if (errcode == c_io_error) return errcode end
        else # Parameter sweep

            # Set parameter lists
            p = Vector[]
            n_dims = output.n_parameters
            for i in 1:n_dims
                if (output.scale[i] == o_scale_lin)
                    push!(p, range(output.x_min[i], output.x_max[i], length=output.x_steps[i]))
                elseif (output.scale[i] == o_scale_log)
                    push!(p, LogRange(output.x_min[i], output.x_max[i], output.x_steps[i]))
                end
            end

            # Initialize step counter
            step = ones(Int64, n_dims)

            # Initialize parameters
            param = zeros(Float64, n_dims)

            # Start loop
            while step[n_dims] <= output.x_steps[n_dims]

                # Update parameter values
                for i in 1:n_dims
                    param[i] = p[i][step[i]]
                end

                # Copy species & reaction lists and system
                species_list_run = copy_species_list(species_list)
                reaction_list_run = copy_reaction_list(reaction_list)
                system_run = copy_system(system)

                errcode = UpdateOutputParameters!(species_list_run, reaction_list_run,
                    system_run, sID, output, param)
                if (errcode == c_io_error) return errcode end

                # Print parameter state
                for i in 1:n_dims
                    @printf("%4s = %10f - ", output.name[i], param[i])
                end
                PrintSpeciesList(species_list_run, sID)

                # Run problem
                sol = @time ExecuteProblem(species_list_run, reaction_list_run, system_run, sID)
                errcode = LoadOutputBlock!(output, sol, species_list_run,
                    reaction_list_run, system_run, sID, param)
                if (errcode == c_io_error) return errcode end

                # Update step
                for i in 1:n_dims
                    step[i] += 1
                    if step[i] > output.x_steps[i] && i<n_dims
                        step[i] = 1
                    else
                        break
                    end
                end

            end # LOOP: over every parameter value
        end # IF: single_run or different
    end # LOOP: Output list

    print("Total time...\n")
    return errcode
end


function UpdateOutputParameters!(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID, output::OutputBlock,
    param::Vector{Float64})
    errcode = 0

    part_press_old = 0.0
    part_press_new = 0.0
    part_press_species = Int64[]
    for i in 1:output.n_parameters
        if (output.case[i] == o_pL)
            # L is kept constant but p is being changed
            # Ideal gas law, p = n kb T, temperature is kept constant
            s_id = output.species_id[i]
            species_list[s_id].pressure = param[i]/system.l
        elseif (output.case[i] == o_dens)
            s_id = output.species_id[i]
            species_list[s_id].pressure = param[i] * kb * species_list[s_id].temp
            # Density is set at the end of this function
        elseif (output.case[i] == o_pressure)
            s_id = output.species_id[i]
            species_list[s_id].pressure = param[i]
        elseif (output.case[i] == o_pressure_percent)
            s_id = output.species_id[i]
            part_press_old += species_list[s_id].pressure / system.total_pressure
            p_partial = param[i] * system.total_pressure 
            species_list[s_id].pressure = p_partial

            # Because partial pressure of this species has changed, the partial
            # pressure of the remaining species need to be readjusted
            part_press_new += param[i]
            push!(part_press_species, s_id)
        elseif (output.case[i] == o_temp)
            s_id = output.species_id[i]
            if s_id == neutral_species_id
                for s in species_list
                    if s.id != sID.electron
                        s.temp = param[i]
                        s.pressure = s.dens * kb * s.temp
                    end
                end
            else
                species_list[s_id].temp = param[i]
            end
        elseif (output.case[i] == o_power)
            system.drivP = param[i]
        end
    end

    # Re-equilibrate pressures
    if (part_press_old != part_press_new)
        part_press_remain_old = 1.0 - part_press_old
        part_press_remain_new = 1.0 - part_press_new
        ratio = part_press_remain_new / part_press_remain_old
        for s in species_list
            account_pressure = true
            # Exclude species whose pressure has been just changed
            for p_id in part_press_species
                if s.id == p_id
                    account_pressure = false
                end
            end

            # Exclude plasma particles
            if s.charge != 0
                account_pressure = false
            end

            if account_pressure
                s.pressure *= ratio
            end
        end
    end

    # Update density value on all species
    for s in species_list
        s.dens = s.pressure / (kb * s.temp)
    end

    return errcode
end


function LoadOutputBlock!(output::OutputBlock, sol,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID, param::Vector{Float64}=Float64[])

    errcode = 0
    print("Load output data...\n")

    n_species = length(species_list)
    if output.case[1] == o_single_run
        output.n_data_frame.time = sol.t
        output.T_data_frame.time = sol.t
        n_steps = length(sol.t)

        # Dump dens/temp into output block
        for s in species_list
            if s.has_dens_eq
                output.n_data_frame[!,s.name] = sol[s.id+n_species, :]
            end
            if s.has_temp_eq
                output.T_data_frame[!,s.name] = sol[s.id, : ]
            end
        end

        # Initialize temp/dens array 
        temp = zeros(n_species)
        dens = zeros(n_species)
        for s in species_list
            if !s.has_temp_eq
                temp[s.id] = s.temp 
            end
            if !s.has_dens_eq
                dens[s.id] = s.dens
            end
        end

        # Loop over every time step
        for i in 1:n_steps
            # Update temperature values
            for s in species_list
                if s.has_temp_eq
                    temp[s.id] = output.T_data_frame[i,s.name]
                end
                if s.has_dens_eq
                    dens[s.id] = output.n_data_frame[i,s.name]
                end
            end

            # Get K values for the curren time step
            K_list = Float64[sol.t[i]]
            for r in reaction_list
                if r.case == r_wall_loss
                    continue
                    #K = r.rate_coefficient(temp, species_list, system, sID) 
                    # Wall loss reactions require to update species parameters
                    # such as h_R, h_L, D, gamma, etc.
                else
                    if system.prerun
                        #push!(K_list, r.rate_coefficient(temp, sID))
                        K = r.rate_coefficient(temp, sID)
                    else
                        #push!(K_list, ReplaceExpressionValues(r.rate_coefficient, temp,
                        #    species_list, system, sID))
                        K = ReplaceExpressionValues(r.rate_coefficient, temp,
                            species_list, system, sID)
                    end
                    push!(K_list, prod(dens[r.reactant_species])*K )
                end
            end

            # Push data into rate coefficient data_frame
            push!(output.K_data_frame, K_list)
        end
    else

        # Initialize buffer lists
        dens_list = copy(param)
        temp_list = copy(param)
        K_list = copy(param)
        temp = zeros(n_species) # temp array is used later for getting rate coefficient values
        dens = zeros(n_species) # dens array is used later for getting rate coefficient values

        # Dump dens/temp into buffer lists 
        for s in species_list
            if s.has_dens_eq
                push!(dens_list, sol[s.id+n_species, end])
            end
            if s.has_temp_eq
                push!(temp_list, sol[s.id, end])
            end
            temp[s.id] = s.temp 
            dens[s.id] = s.dens
        end
        # Push dens/temp buffer lists into data_frames
        push!(output.T_data_frame, temp_list)
        push!(output.n_data_frame, dens_list)

        # Dump K values into buffer 
        for r in reaction_list
            if r.case == r_wall_loss
                continue
                #K = r.rate_coefficient(temp, species_list, system, sID) 
                # Wall loss reactions require to update species parameters
                # such as h_R, h_L, D, gamma, etc.
            else
                if system.prerun
                    #push!(K_list, r.rate_coefficient(temp, sID))
                    K = r.rate_coefficient(temp, sID)
                else
                    #push!(K_list, ReplaceExpressionValues(r.rate_coefficient, temp,
                    #    species_list, system, sID))
                    K = ReplaceExpressionValues(r.rate_coefficient, temp,
                        species_list, system, sID)
                end
                push!(K_list, prod(dens[r.reactant_species])*K )
            end
        end
        # Push K buffer list into rate coefficient data_frame
        push!(output.K_data_frame, K_list)
    end

    return errcode
end


function LogRange(x_min::Float64, x_max::Float64, steps::Int64)
    return [10^y for y in range(log10(x_min), log10(x_max), length=steps)]
end


function copy_system(system::System)

    s = System()

    s.A = system.A
    s.V = system.V
    s.l = system.l
    s.radius = system.radius

    s.power_input_method = system.power_input_method
    s.Vsheath_solving_method = system.Vsheath_solving_method
    s.drivf = system.drivf
    s.drivOmega = system.drivOmega
    s.drivP = system.drivP
    s.P_shape = system.P_shape
    s.P_duty_ratio = system.P_duty_ratio

    s.t_end = system.t_end
    s.alpha = system.alpha
    s.Lambda = system.Lambda

    s.prerun = system.prerun
    s.folder = system.folder
    s.total_pressure = system.total_pressure

    return s
end


function copy_species(s::Species)

    new_s = Species()

    new_s.id = s.id
    new_s.species_id = s.species_id

    new_s.mass = s.mass
    new_s.charge = s.charge

    new_s.has_dens_eq = s.has_dens_eq
    new_s.has_temp_eq = s.has_temp_eq
    new_s.has_wall_loss = s.has_wall_loss
    new_s.has_heating_mechanism = s.has_heating_mechanism
    new_s.has_flow_rate = s.has_flow_rate

    new_s.dens = s.dens
    new_s.temp = s.temp
    new_s.pressure = s.pressure

    new_s.reaction_list = s.reaction_list
    new_s.mfp = s.mfp
    new_s.v_thermal = s.v_thermal
    new_s.v_Bohm = s.v_Bohm
    new_s.D = s.D
    new_s.h_L = s.h_L
    new_s.h_R = s.h_R
    new_s.gamma = s.gamma
    new_s.n_sheath = s.n_sheath
    new_s.flux = s.flux
    new_s.flow_rate = s.flow_rate

    new_s.name = s.name

    return new_s
end


function copy_reaction(r::Reaction)

    new_r = Reaction()

    new_r.name = r.name
    new_r.id = r.id
    new_r.case = r.case
    new_r.neutral_species_id = r.neutral_species_id

    new_r.involved_species = r.involved_species
    new_r.species_balance = r.species_balance
    new_r.reactant_species = r.reactant_species

    new_r.rate_coefficient = r.rate_coefficient

    new_r.E_threshold = r.E_threshold

    return new_r
end


function copy_species_list(species_list::Vector{Species})

    new_species_list = Species[]
    for s in species_list
        new_s = copy_species(s)
        push!(new_species_list, new_s)
    end
    return new_species_list
end


function copy_reaction_list(reaction_list::Vector{Reaction})

    new_reaction_list = Reaction[]
    for r in reaction_list
        new_r = copy_reaction(r)
        push!(new_reaction_list, new_r)
    end
    return new_reaction_list
end

end