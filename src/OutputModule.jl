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
using SharedData: kb
using SharedData: c_io_error, r_emission_rate
using SharedData: o_scale_lin, o_scale_log
using SharedData: o_single_run, o_pL, o_dens, o_temp, o_power, o_pressure
using SharedData: o_pressure_percent, heavy_species_id
using SharedData: o_frequency, o_duty_ratio, o_total_pressure
using EvaluateExpressions: ReplaceExpressionValues
using PlasmaParameters: UpdateParameters!
using PlasmaSheath: GetSheathVoltage!
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using SolveSystem: ExecuteProblem
using PrintModule: PrintErrorMessage, PrintMessage

using CSV
using Printf
using DataFrames: DataFrame
using TimerOutputs

###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GenerateOutputs
#   - SetupInitialConditions
#   - Output_PL

const to = TimerOutput()

function GenerateOutputs!(
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, output_list::Vector{OutputBlock}, sID::SpeciesID)
    
    errcode = 0
    reset_timer!(to)

    for output in output_list
        
        # This flag is used to determine whether headers must be written on
        # new or existing CSV files
        output.first_dump = true

        if output.case[1] == o_single_run
            sol = @timeit to "Problem execution" ExecuteProblem(species_list, reaction_list, system,
                sID, output.case[1])
            errcode = @timeit to "Output dump" LoadOutputBlock!(sol, output, species_list,
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
                system_run = copy_system(system)

                errcode = UpdateOutputParameters!(species_list_run,
                    reaction_list, system_run, sID, output, param)
                if (errcode == c_io_error) return errcode end

                # Print parameter state
                message=""
                for i in 1:n_dims
                    if i == n_dims
                        message *= @sprintf("%s = %g", output.name[i], param[i])
                    else
                        message *= @sprintf("%s = %g; ", output.name[i], param[i])
                    end
                end
                PrintMessage(system_run, message*"\n")

                # Run problem
                sol = @timeit to message ExecuteProblem(species_list_run,
                    reaction_list, system_run, sID)
                errcode = LoadOutputBlock!(sol, output, species_list_run,
                    reaction_list, system_run, sID, param)
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

            PrintMessage(system_run, "\n")
            end # LOOP: over every parameter value
        end # IF: single_run or different
    end # LOOP: Output list

    show(to)
    open(system.log_file,"a") do file
        show(file, to)
    end

    return errcode
end


function UpdateOutputParameters!(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID, output::OutputBlock,
    param::Vector{Float64})
    errcode = 0

    for i in 1:output.n_parameters
        if (output.case[i] == o_pL)
            # L is kept constant but p is being changed
            # Ideal gas law, p = n kb T, temperature is kept constant
            s_id = output.species_id[i]
            s_pressure = param[i]/system.l
            UpdatePressure!(species_list[s_id], s_pressure, system)
        
        elseif (output.case[i] == o_dens)
            s_id = output.species_id[i]
            species_list[s_id].dens = param[i]
            s_pressure = param[i] * kb * species_list[s_id].temp
            UpdatePressure!(species_list[s_id], s_pressure, system)

        elseif (output.case[i] == o_pressure)
            s_id = output.species_id[i]
            s_pressure = param[i]
            UpdatePressure!(species_list[s_id], s_pressure, system)
        
        elseif (output.case[i] == o_total_pressure)
            # Temperature is kept constant but densities are updated
            old_total_press = system.total_pressure
            new_total_press = param[i]
            system.total_pressure = new_total_press 
            for s in species_list
                if s.charge == 0.0
                    p_fract = s.pressure / old_total_press
                    s.pressure = p_fract * new_total_press
                    s.dens = s.pressure / kb / s.temp
                end
            end

        elseif (output.case[i] == o_temp)
            s_id = output.species_id[i]
            if s_id == heavy_species_id
                for s in species_list
                    if s.id != sID.electron
                        s.temp = param[i]
                        #s_pressure = s.dens * kb * s.temp
                        #UpdatePressure!(s, s_pressure, system)
                    end
                end
            else
                species_list[s_id].temp = param[i]
                #s_pressure = species_list[s_id].dens * kb * species_list[s_id].temp
                #UpdatePressure!(species_list[s_id], s_pressure, system)
            end
        elseif (output.case[i] == o_power)
            system.drivP = param[i]
            system.P_absorbed = system.drivP / system.V 
        elseif (output.case[i] == o_frequency)
            system.drivf= param[i]
            system.drivOmega= 2.0*pi*param[i]
        elseif (output.case[i] == o_duty_ratio)
            system.P_duty_ratio = param[i]
        end
    end

    errcode = UpdateOutputParameters_partial_pressure!(species_list, system,
        output, param)
    if errcode != 0
        err_message = "in UpdateOutputParameters while updating partial pressures\n"
        PrintErrorMessage(err_message)
        return c_io_error
    end

    errcode = UpdateOutputParameters_check_pressure(species_list, system)
    if errcode != 0
        err_message = "in UpdateOutputParameters while checking pressures\n"
        PrintErrorMessage(err_message)
        return c_io_error
    end

    return errcode
end


function UpdateOutputParameters_partial_pressure!(species_list::Vector{Species},
    system::System, output::OutputBlock, param::Vector{Float64})

    errcode = 0

    # Change species concentrations (via partial pressures)
    curr_total_pressure = copy(system.total_pressure)
    press_fract_old = 0.0
    press_fract_new = 0.0
    part_fract_species = Int64[]
    for i in 1:output.n_parameters
        if (output.case[i] == o_pressure_percent)
            s_id = output.species_id[i]
            press_fract_old += species_list[s_id].pressure / curr_total_pressure 
            p_partial = param[i] * curr_total_pressure 
            species_list[s_id].pressure = p_partial
            # Because partial pressure of this species has changed, the partial
            # pressure of the remaining species need to be readjusted
            press_fract_new += param[i]
            push!(part_fract_species, s_id)
        end
    end

    # Readjust partial pressures
    if (press_fract_old != press_fract_new)
        press_fract_remain_old = 1.0 - press_fract_old
        press_fract_remain_new = 1.0 - press_fract_new
        ratio = press_fract_remain_new / press_fract_remain_old
        for s in species_list
            press_check = true
            # Exclude species whose pressure fraction has changed
            for p_id in part_fract_species
                if s.id == p_id
                    s.pressure *= system.total_pressure / curr_total_pressure 
                    press_check = false 
                end
            end

            # Exclude plasma particles
            if s.charge != 0
                press_check = false 
            end

            # Species whose pressure is affected because other species
            # fraction has been adjusted
            if press_check
                s.pressure *= ratio
            end
        end
    end

    return errcode
end


function UpdateOutputParameters_check_pressure(species_list::Vector{Species},
    system::System)
    
    errcode = 0

    # Update density value on all species
    # - check total pressure
    press_buffer = 0.0
    for s in species_list
        if isnan(s.pressure)
            err_message = @sprintf("Bad output specification: Species %s pressure is NaN\n", s.name)
            PrintErrorMessage(err_message)
            return c_io_error
        end
        if isinf(s.pressure) 
            err_message = @sprintf("Bad output specification: Species %s pressure is Inf\n", s.name)
            PrintErrorMessage(err_message)
            return c_io_error
        end

        # Update density values
        s.dens = s.pressure / (kb * s.temp)
         if s.charge == 0
            press_buffer += s.pressure
         end
    end

    if abs(press_buffer - system.total_pressure) > 1.e-10
        err_message = "Bad output specification: Sum of species pressure does not match system total pressure\n"
        PrintErrorMessage(err_message)
        return c_io_error
    end

    return errcode
end


function UpdatePressure!(s::Species, s_pressure::Float64, system::System)

    # Updates species pressure and, if not an ion, total_pressure
    if (s.charge == 0)
        system.total_pressure -= s.pressure 
    end
    s.pressure = s_pressure
    if (s.charge == 0)
        system.total_pressure += s.pressure
    end

end


function LoadOutputBlock!(sol, output::OutputBlock,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID, param::Vector{Float64}=Float64[])

    errcode = 0
    PrintMessage(system, "Loading output data...\n")
    T_filename = string(system.folder,"T",output.label,".csv")
    n_filename = string(system.folder,"n",output.label,".csv")
    K_filename = string(system.folder,"K",output.label,".csv")
    param_filename = string(system.folder,"param",output.label,".csv")
    if isfile(T_filename)
        if output.first_dump
            write_header = true
            output.first_dump = false
        else
            write_header = false
        end
    else
        write_header = true
        output.first_dump = false
    end

    n_species = length(species_list)

    # Get neutral/ions temperature
    neutral_temp = -1
    for s in species_list
        if s.id != sID.electron
            neutral_temp = s.temp
            break 
        end
    end

    if output.case[1] == o_single_run
        # - dens & temp data frames are filled by columns:
        #   time and followed by dens/temp species
        # - rate coeff. (K) and parameters data frames are filled by rows:
        #   each time step is a new row
        output.n_data_frame.time = sol.t
        output.T_data_frame.time = sol.t
        n_steps = length(sol.t)

        # Dump dens/temp into output block
        for s in species_list
            output.n_data_frame[!,s.name] = sol[s.id+1,:]
            if s.id == sID.electron 
                output.T_data_frame[!,s.name] = sol[1,:]
            end
        end

        data_len = length(sol[1,:])
        output.T_data_frame[!,"neutrals"] = ones(data_len) * neutral_temp

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
            time = sol.t[i]

            #### Update temperature values
            for s in species_list
                if s.has_temp_eq
                    temp[s.id] = output.T_data_frame[i,s.name]
                end
                if s.has_dens_eq
                    dens[s.id] = output.n_data_frame[i,s.name]
                end
            end

            #### Update plasma parameters
            errcode = UpdateParameters!(temp, dens, species_list, reaction_list, system, sID)
            if errcode == c_io_error
                PrintErrorMessage(system, "UpdateSpeciesParameters  at OutputModule failed")
                return errcode
            end

            #### Update flux and potential values 
            errcode = UpdatePositiveFlux!(species_list)
            if errcode == c_io_error
                PrintErrorMessage(system, "UpdatePositiveFlux failed")
                return errcode
            end
            errcode = GetSheathVoltage!(system, species_list, sID, time)
            if errcode == c_io_error
                PrintErrorMessage(system, "GetSheathVoltage failed")
                return errcode
            end
            errcode = UpdateNegativeFlux!(species_list, system, sID)
            if errcode == c_io_error
                PrintErrorMessage(system, "UpdateNegativeFlux failed")
                return errcode
            end
            
            # Get K values for the curren time step
            K_list = Float64[time]
            for r in reaction_list
                push!(K_list, prod(dens[r.reactant_species])*r.K_value )
            end
            # Push K data into rate coefficient data_frame
            push!(output.K_data_frame, K_list)

            # Gather and push parameter data into parameter data_frame
            param_list = Float64[time, system.plasma_potential]
            # Push total pressure
            push!(param_list, system.total_pressure)
            # Push electronegativity 
            push!(param_list, system.alpha)

            for s in species_list
                if !(s.charge==0)
                    push!(param_list, s.flux)
                    push!(param_list, s.mfp)
                end
            end
            push!(output.param_data_frame, param_list)

        end
        CSV.write(T_filename, output.T_data_frame)
        CSV.write(n_filename, output.n_data_frame)
        CSV.write(K_filename, output.K_data_frame)
        CSV.write(param_filename, output.param_data_frame)
    else

        # Initialize buffer lists
        dens_list = copy(param)
        temp_list = copy(param)
        K_list = copy(param)
        p_list = copy(param)
        temp = zeros(n_species) # temp array is used later for getting rate coefficient values
        dens = zeros(n_species) # dens array is used later for getting rate coefficient values

        # Various parameter dump
        push!(p_list, system.plasma_potential)
        push!(p_list, system.total_pressure)
        push!(p_list, system.alpha)

        # Dump dens/temp and other into buffer lists 
        for s in species_list
            push!(dens_list, sol[s.id+1, end])
            if s.id == sID.electron
                push!(temp_list, sol[1, end])
            end
            temp[s.id] = s.temp 
            dens[s.id] = s.dens

            if !(s.charge==0)
                push!(p_list, s.flux)
                push!(p_list, s.mfp)
            end
        end
        # Add neutral temperature to temperature array
        push!(temp_list, neutral_temp)

        # Push dens/temp buffer lists into data_frames
        push!(output.T_data_frame, temp_list)
        push!(output.n_data_frame, dens_list)
        push!(output.param_data_frame, p_list)
        CSV.write(T_filename, DataFrame(output.T_data_frame[end,:]),
            append=true, writeheader=write_header)
        CSV.write(n_filename, DataFrame(output.n_data_frame[end,:]),
            append=true, writeheader=write_header)
        CSV.write(param_filename, DataFrame(output.param_data_frame[end,:]),
            append=true, writeheader=write_header)

        # Dump K values into buffer 
        for r in reaction_list
            push!(K_list, prod(dens[r.reactant_species])*r.K_value )
        end
        # Push K buffer list into rate coefficient data_frame
        push!(output.K_data_frame, K_list)
        CSV.write(K_filename, DataFrame(output.K_data_frame[end,:]),
            append=true , writeheader=write_header)
    end

    PrintMessage(system, "Loading data done\n")
    return errcode
end


function LogRange(x_min::Float64, x_max::Float64, steps::Int64)
    return [10^y for y in range(log10(x_min), log10(x_max), length=steps)]
end


function copy_system(system::System)

    s = System()

    s.A = copy(system.A)
    s.V = copy(system.V)
    s.l = copy(system.l)
    s.radius = copy(system.radius)

    s.h_id = copy(system.h_id)
    s.Vsheath_solving_method = copy(system.Vsheath_solving_method)
    s.drivf = copy(system.drivf)
    s.drivOmega = copy(system.drivOmega)
    s.drivP = copy(system.drivP)
    s.P_absorbed = copy(system.P_absorbed)
    s.P_shape = copy(system.P_shape)
    s.P_duty_ratio = copy(system.P_duty_ratio)
    s.P_start = copy(system.P_start)
    s.dt_start = copy(system.dt_start)

    s.plasma_potential = copy(system.plasma_potential)
    s.total_pressure = copy(system.total_pressure)
    s.T_e_min = copy(system.T_e_min)
    s.T_e_max = copy(system.T_e_max)

    s.t_end = copy(system.t_end)
    s.errcode = copy(system.errcode)

    s.alpha = copy(system.alpha)
    s.Lambda = copy(system.Lambda)

    s.prerun = copy(system.prerun)
    s.folder = system.folder
    s.log_file = system.log_file

    return s
end


function copy_species(s::Species)

    new_s = Species()

    new_s.id = copy(s.id)
    new_s.species_id = copy(s.species_id)

    new_s.mass = copy(s.mass)
    new_s.charge = copy(s.charge)

    new_s.has_dens_eq = copy(s.has_dens_eq)
    new_s.has_temp_eq = copy(s.has_temp_eq)
    new_s.has_wall_loss = copy(s.has_wall_loss)
    new_s.has_heating_mechanism = copy(s.has_heating_mechanism)
    new_s.has_flow_rate = copy(s.has_flow_rate)

    new_s.dens = copy(s.dens)
    new_s.temp = copy(s.temp)
    new_s.pressure = copy(s.pressure)

    new_s.reaction_list = s.reaction_list
    new_s.mfp = copy(s.mfp)
    new_s.v_thermal = copy(s.v_thermal)
    new_s.v_Bohm = copy(s.v_Bohm)
    new_s.n_sheath = copy(s.n_sheath)
    new_s.flux = copy(s.flux)
    new_s.D = copy(s.D)
    new_s.h_R = copy(s.h_R)
    new_s.h_L = copy(s.h_L)
    new_s.gamma = copy(s.gamma)
    new_s.flow_rate = copy(s.flow_rate)

    new_s.in_nodes = Tuple{Int64, String, Float64}[]
    for in_tuple in s.in_nodes
        push!(new_s.in_nodes, (in_tuple[1], in_tuple[2], in_tuple[3]))
    end
    
    new_s.out_nodes = Tuple{Int64, String, Float64}[]
    for out_tuple in s.out_nodes
        push!(new_s.out_nodes, (out_tuple[1], out_tuple[2], out_tuple[3]))
    end

    new_s.name = s.name

    return new_s
end


function copy_reaction(r::Reaction)

    new_r = Reaction()

    new_r.name = r.name
    new_r.id = copy(r.id)
    new_r.case = copy(r.case)
    new_r.neutral_species_id = copy(r.neutral_species_id)

    new_r.involved_species = copy(r.involved_species)
    new_r.species_balance = copy(r.species_balance)
    new_r.reactant_species = copy(r.reactant_species)
    new_r.product_species = copy(r.product_species)

    new_r.rate_coefficient = r.rate_coefficient
    new_r.K_value = copy(r.K_value)

    new_r.E_threshold = copy(r.E_threshold)

    new_r.self_absorption = copy(r.self_absorption)
    new_r.g_high = copy(r.g_high)
    new_r.g_low = copy(r.g_low)
    new_r.g_high_total = copy(r.g_high_total)
    new_r.g_low_total = copy(r.g_low_total)
    new_r.wavelength = copy(r.wavelength)

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