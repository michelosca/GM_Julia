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

        if output.case == o_single_run
            sol = ExecuteProblem(species_list, reaction_list, system, sID)
            errcode = LoadOutputBlock!(output, sol, species_list,
                reaction_list, system, sID)
            if (errcode == c_io_error) return errcode end
        else
            # Setup the parameter list to be looped through
            if (output.scale == o_scale_lin)
                x_array = range(output.x_min, output.x_max, length=output.x_steps)
            elseif (output.scale == o_scale_log)
                x_array = LogRange(output.x_min, output.x_max, output.x_steps)
            end
            
            # Start parameter loop
            for x in x_array
                errcode = UpdateOutputParameters!(species_list, reaction_list,
                    system, sID, output, x)
                if (errcode == c_io_error) return errcode end

                @printf("%4s = %10f - ", output.parameter, x)
                sol = ExecuteProblem(species_list, reaction_list, system, sID)
                errcode = LoadOutputBlock!(output, sol, species_list,
                    reaction_list, system, sID, x)
                if (errcode == c_io_error) return errcode end
            end
        end
    end # Output list loop

    return errcode
end


function UpdateOutputParameters!(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID, output::OutputBlock,
    x::Float64)
    errcode = 0

    if (output.case == o_pL)
        # L is kept constant but p is being changed
        # Ideal gas law, p = n kb T, temperature is kept constant
        p = x / system.l 
        s_id = output.species_id
        T = species_list[s_id].temp
        species_list[s_id].dens = p / kb / T
    elseif (output.case == o_dens)
        s_id = output.species_id
        species_list[s_id].dens = x 
    elseif (output.case == o_pressure)
        s_id = output.species_id
        T = species_list[s_id].temp
        species_list[s_id].dens = x / kb / T
    elseif (output.case == o_temp)
        s_id = output.species_id
        species_list[s_id].temp = x 
    elseif (output.case == o_power)
        system.drivP = x
    end



    return errcode
end


function LoadOutputBlock!(output::OutputBlock, sol,
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID, param::Float64 = 0.0)

    errcode = 0

    n_species = length(species_list)
    if output.case == o_single_run
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
            # Update temperature values for each species
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
            push!(output.K_data_frame, K_list)
        end
    else

        # Initialize buffer lists
        dens_list = Float64[param]
        temp_list = Float64[param]
        K_list = Float64[param]
        temp = zeros(n_species)

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
        push!(output.K_data_frame, K_list)
    end

    return errcode
end

function LogRange(x_min::Float64, x_max::Float64, steps::Int64)
    return [10^y for y in range(log10(x_min), log10(x_max), length=steps)]
end

end