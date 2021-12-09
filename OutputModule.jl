module OutputModule

using SharedData: Species, Reaction, System, SpeciesID, OutputBlock
using SharedData: kb, K_to_eV
using SharedData: c_io_error, r_wall_loss
using SharedData: o_scale_lin, o_scale_log
using SharedData: o_single_run, o_pL, o_dens, o_temp, o_power, o_pressure
using EvaluateExpressions: ReplaceExpressionValues
using SolveSystem: ExecuteProblem
using Printf

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
            sol = ExecuteProblem(species_list, reaction_list, system, sID)
            errcode = LoadOutputBlock!(output, sol, species_list,
                reaction_list, system, sID)
            if (errcode == c_io_error) return errcode end
        else
            n_dims = output.n_parameters

            # Set parameter lists
            p = Vector[]
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
                errcode = UpdateOutputParameters!(species_list, reaction_list,
                    system, sID, output, param)
                if (errcode == c_io_error) return errcode end

                # Print parameter state
                for i in 1:n_dims
                    @printf("%4s = %10f - ", output.name[i], param[i])
                end

                # Run problem
                sol = ExecuteProblem(species_list, reaction_list, system, sID)
                errcode = LoadOutputBlock!(output, sol, species_list,
                    reaction_list, system, sID, param)
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
            p = param[i] / system.l 
            s_id = output.species_id[i]
            T = species_list[s_id].temp
            species_list[s_id].dens = p / kb / T
        elseif (output.case[i] == o_dens)
            s_id = output.species_id[i]
            species_list[s_id].dens = param[i]
        elseif (output.case[i] == o_pressure)
            s_id = output.species_id[i]
            T = species_list[s_id].temp
            species_list[s_id].dens = param[i] / kb / T
        elseif (output.case[i] == o_temp)
            s_id = output.species_id[i]
            species_list[s_id].temp = param[i]
        elseif (output.case[i] == o_power)
            system.drivP = param[i]
        end
    end

    return errcode
end


function LoadOutputBlock!(output::OutputBlock, sol,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID, param::Vector{Float64}=Float64[])

    errcode = 0

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

        # Initialize temp array 
        temp = zeros(n_species)
        for s in species_list
            if !s.has_temp_eq
                temp[s.id] = s.temp 
            end
        end

        # Loop over every time step
        for i in 1:n_steps
            # Update temperature values
            for s in species_list
                if s.has_temp_eq
                    temp[s.id] = output.T_data_frame[i,s.name]
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
                        push!(K_list, r.rate_coefficient(temp, sID))
                    else
                        push!(K_list, ReplaceExpressionValues(r.rate_coefficient, temp,
                            species_list, system, sID))
                    end
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

        # Dump dens/temp into buffer lists 
        for s in species_list
            if s.has_dens_eq
                push!(dens_list, sol[s.id+n_species, end])
            end
            if s.has_temp_eq
                push!(temp_list, sol[s.id, end])
            end
            temp[s.id] = s.temp 
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
                    push!(K_list, r.rate_coefficient(temp, sID))
                else
                    push!(K_list, ReplaceExpressionValues(r.rate_coefficient, temp,
                        species_list, system, sID))
                end
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

end