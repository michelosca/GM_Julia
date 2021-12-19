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
                if output.x_max[i] <= output.x_min[i]
                    print("***ERROR*** Output block must fulfil that x_max > x_max\n")
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
                end
            elseif (output.case[i] == o_dens)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for dens outputs\n")
                end
            elseif (output.case[i] == o_pressure)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for pressure outputs\n")
                end
            elseif (output.case[i] == o_temp)
                if output.species_id[i] == 0
                    print("***ERROR*** Must specify the 'species' entry for temp outputs\n")
                end
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
    output = output_list[end]

    if (name == "parameter" || name == "x")
        SetUpEmptyParameter!(output)
        i = output.n_parameters
        if (var=="pL")
            output.case[i] = o_pL
        elseif (var=="dens" || var=="density")
            output.case = o_dens
        elseif (var=="temp" || var=="temperature")
            output.case[i] = o_temp
        elseif (var=="pressure")
            output.case[i] = o_pressure
        elseif (var=="input_power" || var=="power")
            output.case[i] = o_power
        elseif (var=="time" || var=="t")
            output.case[i] = o_single_run
        end
    end


    i = output.n_parameters
    if (i == 0)
        print("***ERROR*** The first declaration in output block must be 'parameter' or 'x'\n")
        return c_io_error
    end

    if (name == "species")
        if (i < length(output.species_id))
            print("***ERROR*** Output block: species has been declared before 'parameter'/'x'\n")
            return c_io_error
        end
        for s in species_list
            if var == s.name
                output.species_id[i] = s.id
                break
            end
        end
    end
    
    if (name == "x_min")
        if (i < length(output.x_min))
            print("***ERROR*** Output block: x_min has been declared before 'parameter'/'x'\n")
            return c_io_error
        end
        output.x_min[i] = parse(Float64, var) * units
    end

    if (name == "x_max")
        if (i < length(output.x_max))
            print("***ERROR*** Output block: x_max has been declared before 'parameter'/'x'\n")
            return c_io_error
        end
        output.x_max[i] = parse(Float64, var) * units
    end

    if (name == "steps")
        if (i < length(output.x_steps))
            print("***ERROR*** Output block: x_steps has been declared before 'parameter'/'x'\n")
            return c_io_error
        end
        output.x_steps[i] = parse(Int64, var)
    end

    if (name == "scale")
        if (i < length(output.scale))
            print("***ERROR*** Output block: scale has been declared before 'parameter'/'x'\n")
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
    output.name = String[]
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
            sname = species_list[output.species_id[i]].name
            output.name[i] = string("T_",sname)
        elseif output.case[i] == o_pressure
            sname = species_list[output.species_id[i]].name
            output.name[i] = string("P_",sname)
        elseif output.case[i] == o_single_run
            output.name[i] = "time"
        else
            errcode = c_io_error
            print("***ERROR*** Output parameter not recognized\n")
        end

        if output.case[i] != o_single_run
            output.n_data_frame[!, output.name[i]] = Float64[]
            output.T_data_frame[!, output.name[i]] = Float64[]
            if i == output.n_parameters 
                for s in species_list
                    if s.has_dens_eq
                        output.n_data_frame[!, s.name] = Float64[]
                    end
                    if s.has_temp_eq
                        output.T_data_frame[!, s.name] = Float64[]
                    end
                end
            end
        end

        # Initialize rate coefficient K_data_frame
        output.K_data_frame[!, output.name[i]] = Float64[]
        if i == output.n_parameters
            for r in reaction_list
                if r.case == r_wall_loss
                    continue
                end
                col_label = string("r",r.id,": ",r.name)
                output.K_data_frame[!, col_label] = Float64[]
            end
        end
    end

    return errcode
end

end