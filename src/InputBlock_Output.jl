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

module InputBlock_Output

using SharedData: OutputBlock, Species, Reaction, System
using SharedData: c_io_error
using SharedData: o_scale_lin, o_scale_log
using SharedData: r_wall_loss
using SharedData: o_single_run, o_pL, o_dens, o_temp, o_power, o_pressure
using SharedData: o_pressure_percent, neutral_species_id
using SharedData: o_frequency, o_duty_ratio, o_total_pressure
using InputBlock_System: GetUnits!
using DataFrames: DataFrame
using PrintModule: PrintErrorMessage

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
            SetUpEmptyParameter!(output)
            output.case[1] = o_single_run
    
            errcode = SetupOutputBlock!(output, species_list, reaction_list)

            # Add output block to output_list
            push!(output_list, output)
        else
            for output in output_list
                errcode = SetupOutputBlock!(output, species_list, reaction_list)
            end
        end

        # Setup output labels
        for output in output_list
            for i in 1:output.n_parameters
                output.label = string(output.label,"_vs_",output.name[i])
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

function EndOutputBlock!(read_step::Int64, output_list::Vector{OutputBlock},
    system::System)
    
    errcode = 0

    if (read_step == 2)
        output = output_list[end]
        n_dims = output.n_parameters
        for i in 1:n_dims
            if output.case[i] == 0
                print("***ERROR*** Output block: Parameter id(case) has not been defined\n")
                errcode = c_io_error
            end

            if output.case[i] == o_single_run && n_dums > 1
                print("***ERROR*** Output block: Single run output does not make sense with several parameter declarations\n")
                errcode = c_io_error
            end

            if output.case[i] != o_single_run
                if output.x_max[i] <= output.x_min[i] && output.x_steps[i] > 1
                    print("***ERROR*** Output block must fulfil that x_max > x_mix\n")
                    errcode = c_io_error
                end

                if output.x_steps[i] <= 0
                    print("***ERROR*** Output block parameter steps <= 0\n")
                    errcode = c_io_error
                end
            end

            if (output.case[i] == o_pL)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for pL outputs\n")
                    errcode = c_io_error
                end
            end
            if (output.case[i] == o_dens)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for dens outputs\n")
                    errcode = c_io_error
                end
            end
            if (output.case[i] == o_pressure || output.case[i] == o_pressure_percent)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for pressure outputs\n")
                    errcode = c_io_error
                end
            end

            if output.case[i] == o_pressure_percent
                if (system.total_pressure <= 0)
                    print("***ERROR*** Partial pressure output requires an in advance 'total_pressure' declaration in SYSTEM block\n")
                    errcode = c_io_error
                end
            end

            if (output.case[i] == o_temp)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for temp outputs\n")
                    errcode = c_io_error
                end
            end

            if (output.case[i] == o_frequency)
                if system.P_shape != "square"
                    print("***ERROR*** P_shape must be 'square' in order to use 'frequency' as output parameter\n")
                    errcode = c_io_error
                end
                if system.P_duty_ratio == 1.0
                    print("***WARNING*** duty_ratio is set to 1.0\n")
                end
            end

            if (output.case[i] == o_duty_ratio)
                if system.P_shape != "square"
                    print("***ERROR*** P_shape must be 'square' in order to use 'duty_ratio' as output parameter\n")
                    errcode = c_io_error
                end

                if output.x_min[i] <= 0.0
                    print("***ERROR*** Minimum duty ratio must be non-zero and positive\n")
                    errcode = c_io_error
                end
                if output.x_max[i] > 1.0
                    print("***ERROR*** Maximum duty ratio must be lower or equal 1\n")
                    errcode = c_io_error
                end
            end
        end
    end

    return errcode
end

function ReadOutputEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, output_list::Vector{OutputBlock},
    species_list::Vector{Species}, system::System)
    
    errcode = 0

    if (read_step == 1)
        return errcode
    end

    units, name = GetUnits!(name)
    output = output_list[end]

    if (name == "parameter" || name == "x")
        SetUpEmptyParameter!(output)
        i = output.n_parameters
        if (var=="pL")
            output.case[i] = o_pL
        elseif (var=="dens" || var=="density")
            output.case[i] = o_dens
        elseif (var=="temp" || var=="temperature")
            output.case[i] = o_temp
        elseif (var=="pressure")
            output.case[i] = o_pressure
        elseif (var=="pressure_%" || var=="pressure_percent" || var=="partial_pressure")
            output.case[i] = o_pressure_percent
        elseif (var=="input_power" || var=="power")
            output.case[i] = o_power
        elseif (var=="time" || var=="t")
            output.case[i] = o_single_run
        elseif (var=="frequency" || var=="f")
            output.case[i] = o_frequency
        elseif (var=="duty_ratio" || var=="dutyratio")
            output.case[i] = o_duty_ratio
        elseif (var=="total_pressure" || var=="p_total")
            output.case[i] = o_total_pressure
        else
            PrintErrorMessage(system, "Output parameter was not found")
            return c_io_error
        end
    end


    i = output.n_parameters
    if (i == 0)
        PrintErrorMessage(system, "The first declaration in output block must be 'parameter' or 'x'")
        return c_io_error
    end

    if (name == "species")
        if (i < length(output.species_id))
            PrintErrorMessage(system, "Output block: species has been declared before 'parameter'/'x'")
            return c_io_error
        end
        for s in species_list
            if var == s.name
                output.species_id[i] = s.id
                break
            end
        end
        if var == "neutral_species" || var == "neutral_gas"
            output.species_id[i] = neutral_species_id
        end
    end
    
    if (name == "x_min")
        if (i < length(output.x_min))
            PrintErrorMessage(system, "Output block: x_min has been declared before 'parameter'/'x'")
            return c_io_error
        end
        output.x_min[i] = parse(Float64, var) * units
    end

    if (name == "x_max")
        if (i < length(output.x_max))
            PrintErrorMessage(system, "Output block: x_max has been declared before 'parameter'/'x'")
            return c_io_error
        end
        output.x_max[i] = parse(Float64, var) * units
    end

    if (name == "steps")
        if (i < length(output.x_steps))
            PrintErrorMessage(system, "Output block: x_steps has been declared before 'parameter'/'x'")
            return c_io_error
        end
        output.x_steps[i] = parse(Int64, var)
    end

    if (name == "scale")
        if (i < length(output.scale))
            PrintErrorMessage(system, "Output block: scale has been declared before 'parameter'/'x'")
            return c_io_error
        end
        if (var == "log" || var == "logarithmic")
            output.scale[i] = o_scale_log
        end
    end

    return errcode
end


function InitializeOutputBlock!(output::OutputBlock)
    output.n_parameters = 0
    output.case = Int64[]
    output.species_id = Int64[]
    output.scale = Int64[]
    output.x = Float64[]
    output.x_min = Float64[]
    output.x_max = Float64[]
    output.x_steps = Int64[] 
    output.n_data_frame = DataFrame()
    output.T_data_frame = DataFrame()
    output.K_data_frame = DataFrame()
    output.param_data_frame = DataFrame()
    output.name = String[]
    output.label = ""
    output.first_dump = true
end


function SetUpEmptyParameter!(output::OutputBlock)

    output.n_parameters += 1
    push!(output.case, 0)
    push!(output.species_id, 0)
    push!(output.scale, o_scale_lin)
    push!(output.x, 0.0)
    push!(output.x_min, 0.0)
    push!(output.x_max, 0.0)
    push!(output.x_steps, 0)
    push!(output.name, "empty")

end


function SetupOutputBlock!(output::OutputBlock,
    species_list::Vector{Species}, reaction_list::Vector{Reaction})

    errcode = 0

    for i in 1:output.n_parameters
        if output.case[i] == o_pL
            sname = species_list[output.species_id[i]].name
            output.name[i] = string("pL_",sname)
        elseif output.case[i] == o_power
            output.name[i] = "power"
        elseif output.case[i] == o_dens
            sname = species_list[output.species_id[i]].name
            output.name[i] = string("n_",sname)
        elseif output.case[i] == o_temp
            s_id = output.species_id[i]
            if s_id == neutral_species_id
                sname = "neutrals"
            else
                sname = species_list[s_id].name
            end
            output.name[i] = string("T_",sname)
        elseif output.case[i] == o_pressure
            sname = species_list[output.species_id[i]].name
            output.name[i] = string("P_",sname)
        elseif output.case[i] == o_pressure_percent
            sname = species_list[output.species_id[i]].name
            output.name[i] = string("P%_",sname)
        elseif output.case[i] == o_single_run
            output.name[i] = "time"
        elseif output.case[i] == o_frequency
            output.name[i] = "frequency" 
        elseif output.case[i] == o_duty_ratio
            output.name[i] = "duty_ratio" 
        elseif output.case[i] == o_total_pressure
            output.name[i] = "total_pressure" 
        else
            errcode = c_io_error
            print("***ERROR*** Output parameter not recognized\n")
        end

        if output.case[i] != o_single_run
            output.n_data_frame[!, output.name[i]] = Float64[]
            output.T_data_frame[!, output.name[i]] = Float64[]
            if i == output.n_parameters 
                for s in species_list
                    output.n_data_frame[!, s.name] = Float64[]
                    if s.has_temp_eq
                        output.T_data_frame[!, s.name] = Float64[]
                    end
                end
            end
        end

        output.K_data_frame[!, output.name[i]] = Float64[]
        output.param_data_frame[!, output.name[i]] = Float64[]
        if i == output.n_parameters
            # Initialize rate coefficient K_data_frame
            for r in reaction_list
                col_label = string("r",r.id,": ",r.name)
                output.K_data_frame[!, col_label] = Float64[]
            end

            # Initialize parameter data frame
            output.param_data_frame[!, "V_plasma"] = Float64[]
            output.param_data_frame[!, "P_total"] = Float64[]
            output.param_data_frame[!, "Electronegativity"] = Float64[]
            for s in species_list
                if !(s.charge == 0)
                    flux_label = string(s.name * "_flux")
                    mfp_label = string(s.name * "_mfp")
                    output.param_data_frame[!, flux_label] = Float64[]
                    output.param_data_frame[!, mfp_label] = Float64[]
                end
            end
        end
    end

    return errcode
end

end