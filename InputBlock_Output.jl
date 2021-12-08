module InputBlock_Output

using SharedData: OutputBlock, Species, Reaction
using SharedData: o_scale_lin, o_scale_log
using SharedData: r_wall_loss
using SharedData: o_single_run, o_pL, o_dens, o_temp, o_power, o_pressure
using InputBlock_System: GetUnits!
using DataFrames: DataFrame

function StartFile_Output!(read_step::Int64, output_list::Vector{OutputBlock})
    
    errcode = 0

    return errcode
end

function EndFile_Output!(read_step::Int64, output_list::Vector{OutputBlock},
    species_list::Vector{Species}, reaction_list::Vector{Reaction})
    
    errcode = 0

    if (read_step == 2)
        if output_list == OutputBlock[]
            # In case no output block is defined, create one for
            # GM time dependent results
            output = OutputBlock()
            InitializeOutputBlock!(output)
            output.case = o_single_run
    
            errcode = InitializeOutputBlockVectors!(output, species_list, reaction_list)

            # Add output block to output_list
            push!(output_list, output)
        else
            for output in output_list
                errcode = InitializeOutputBlockVectors!(output, species_list, reaction_list)
            end
        end
    end
    return errcode
end

function StartOutputBlock!(read_step::Int64, output_list::Vector{OutputBlock})
    
    errcode = 0

    if (read_step == 2)
        output = OutputBlock()
        InitializeOutputBlock!(output)
        push!(output_list, output)
    end

    return errcode
end

function EndOutputBlock!(read_step::Int64, output_list::Vector{OutputBlock})
    
    errcode = 0

    if (read_step == 2)
        current_output = output_list[end]

        if current_output.case != o_single_run
            if current_output.x_max <= current_output.x_min
                print("***ERROR*** Output block must fulfil that x_max > x_max\n")
                errcode = c_io_error
            end

            if current_output.x_steps <= 0
                print("***ERROR*** Output block parameter steps <= 0\n")
                errcode = c_io_error
            end

        end

        if (current_output.case == o_pL)
            if current_output.species_id == 0
                print("***ERROR*** Must specify the 'species' entry for pL outputs\n")
            end
        elseif (current_output.case == o_dens)
            if current_output.species_id == 0
                print("***ERROR*** Must specify the 'species' entry for dens outputs\n")
            end
        elseif (current_output.case == o_pressure)
            if current_output.species_id == 0
                print("***ERROR*** Must specify the 'species' entry for pressure outputs\n")
            end
        elseif (current_output.case == o_temp)
            if current_output.species_id == 0
                print("***ERROR*** Must specify the 'species' entry for temp outputs\n")
            end
        end

    end

    return errcode
end

function ReadOutputEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, output_list::Vector{OutputBlock},
    species_list::Vector{Species})
    
    errcode = 0

    if (read_step == 1)
        return errcode
    end

    units, name = GetUnits!(name)

    if (name == "parameter" || name == "x")
        if (var=="pL")
            output_list[end].case = o_pL
        elseif (var=="dens" || var=="density")
            output_list[end].case = o_dens
        elseif (var=="temp" || var=="temperature")
            output_list[end].case = o_temp
        elseif (var=="pressure")
            output_list[end].case = o_pressure
        elseif (var=="input_power" || var=="power")
            output_list[end].case = o_power
        end
    end

    if (name == "species")
        for s in species_list
            if var == s.name
                output_list[end].species_id = s.id
                break
            end
        end
    end
    
    if (name == "x_min")
        output_list[end].x_min = parse(Float64, var) * units
    end

    if (name == "x_max")
        output_list[end].x_max = parse(Float64, var) * units
    end

    if (name == "steps")
        output_list[end].x_steps = parse(Int64, var)
    end

    if (name == "scale")
        if (var == "log" || var == "logarithmic")
            output_list[end].scale = o_scale_log
        end
    end

    return errcode
end


function InitializeOutputBlock!(outputblock::OutputBlock)
    outputblock.case = 0
    outputblock.species_id = 0
    outputblock.scale = o_scale_lin
    outputblock.x = 0.0 
    outputblock.x_min = 0.0
    outputblock.x_max = 0.0
    outputblock.x_steps = 0
    outputblock.n_data_frame = DataFrame()
    outputblock.T_data_frame = DataFrame()
    outputblock.K_data_frame = DataFrame()
end


function InitializeOutputBlockVectors!(output::OutputBlock,
    species_list::Vector{Species}, reaction_list::Vector{Reaction})

    errcode = 0

    if output.case == o_pL
        output.parameter = "pL"
    elseif output.case == o_power
        output.parameter = "P"
    elseif output.case == o_dens
        output.parameter = "n"
    elseif output.case == o_temp
        output.parameter = "T"
    elseif output.case == o_pressure
        output.parameter = "P"
    elseif output.case == o_single_run
        output.parameter = "time"
    else
        errcode = c_io_error
        print("***ERROR*** Output parameter not recognized\n")
    end

    if output.case != o_single_run
        output.n_data_frame[!, output.parameter] = Float64[]
        output.T_data_frame[!, output.parameter] = Float64[]
        for s in species_list
            if s.has_dens_eq
                output.n_data_frame[!, s.name] = Float64[]
            end
            if s.has_temp_eq
                output.T_data_frame[!, s.name] = Float64[]
            end
        end
    end

    # Initialize rate coefficient K_data_frame
    output.K_data_frame[!, output.parameter] = Float64[]
    for r in reaction_list
        if r.case == r_wall_loss
            continue
        end
        output.K_data_frame[!,string("r",r.id)] = Float64[]
    end

    return errcode
end

end