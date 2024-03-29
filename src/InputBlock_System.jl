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

module InputBlock_System

using SharedData: c_io_error
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation
using SharedData: h_classical, h_Gudmundsson, h_Monahan, h_Thorsteinsson 
using SharedData: p_constant, p_square
using SharedData: System
using SharedData: K_to_eV, e
using PlasmaParameters: GetLambda
using Dates


###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - StartFile_System
# - StartSystemBlock
# - ReadSystemEntry
# - EndSystemBlock
#   - AddSystemToList

function StartFile_System!(read_step::Int64, system::System, filename::String) 

    errcode = 0

    if read_step == 1
        # Generate a log file name
        now_stamp = Dates.now()
        log_num = Dates.format(now_stamp, "yyyymmddHHMMSS")
        log_file = "GM_run_" * log_num * ".log"

        index = findlast("/", filename)
        if index === nothing
            system.folder = "./"
        else
            system.folder = filename[1:index[1]]
        end
        system.log_file = system.folder * log_file
    end

    return errcode
end


function StartSystemBlock!(read_step::Int64, system::System)

    errcode = 0
    if (read_step == 1)
        system.A = 0.0
        system.V = 0.0
        system.l = 0.0
        system.radius = 0.0
        system.h_id = 0
        system.Vsheath_solving_method = 0
        system.drivf = 0.0
        system.drivOmega = 0.0
        system.drivP = 0.0
        system.P_absorbed = 0.0
        system.P_shape = p_constant
        system.P_duty_ratio = 1.0
        system.P_start = 0.0
        system.dt_start = 1.e-12
        system.on_slope = 0.0
        system.off_slope = 0.0
        system.P_off = 0.0
        system.t_end = 0.0
        system.total_pressure = 0.0
        system.T_e_min = 1.5 / K_to_eV
        system.T_e_max = 4.0 / K_to_eV
        system.Lambda = 0.0
        system.K_recombination = 0.0
        system.positive_ion_dens = 0.0
        system.negative_ion_dens = 0.0
        system.prerun = true
        system.plasma_potential = 100.0
        system.errcode = 0
        system.alpha = 0.0
    end
    return errcode
end


function ReadSystemEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, system::System)

    errcode = c_io_error 

    if (read_step == 1)
        
        units_fact, name = GetUnits!(name)

        # Identify variable
        lname = lowercase(name)
        if (name=="A" || lname=="area")
            system.A = parse(Float64,var)
            errcode = 0
        elseif (lname=="v" || lname=="volume" || lname=="vol")
            system.V = parse(Float64, var)
            errcode = 0
        elseif (lname=="l" || lname=="length")
            system.l = parse(Float64, var) 
            errcode = 0
        elseif (lname=="radius" || lname=="r")
            system.radius = parse(Float64, var) 
            errcode = 0
        elseif (name=="f" || lname=="frequency" ||
            lname=="freq" || lname=="driving_frequency")
            system.drivf = parse(Float64,var) * units_fact
            errcode = 0
        elseif (name=="P" || lname=="power_input" ||
            lname=="power")
            system.drivP = parse(Float64, var) * units_fact
            errcode = 0
        elseif (name=="P_shape" || lname=="power_shape")
            if var=="constant"
                system.P_shape = p_constant
                errcode = 0
            elseif var=="square"
                system.P_shape = p_square
                errcode = 0
            else
                print("***ERROR*** Power waveform is not defined\n")
                errcode = c_io_error
            end
        elseif (name=="P_duty_ratio" || lname=="power_duty_ratio")
            system.P_duty_ratio = parse(Float64, var)
            errcode = 0
        elseif (name=="P_t_start" || lname=="power_t_start" ||
            name=="P_time_start" || lname=="power_time_start")
            system.P_start = parse(Float64, var)
            errcode = 0
        elseif (name=="P_off" || lname=="power_off") 
            system.P_off = parse(Float64, var)
            errcode = 0
        elseif (lname=="dt_start")
            system.dt_start = parse(Float64, var)
            errcode = 0
        elseif (lname=="on_slope") || (lname == "p_on_slope")
            system.on_slope = parse(Float64, var)
            errcode = 0
        elseif (lname=="off_slope") || (lname == "p_off_slope")
            system.off_slope = parse(Float64, var)
            errcode = 0
        elseif (lname=="h_factor" || lname=="h" || lname=="h_id")
            lvar = lowercase(var)
            if (lvar == "gudmundsson")
                h_id = h_Gudmundsson
                errcode = 0
            elseif (lvar == "classical")
                h_id = h_classical
                errcode = 0 
            elseif (lvar == "monahan")
                h_id = h_Monahan
                errcode = 0 
            elseif (lvar == "thorsteinsson")
                h_id = h_Thorsteinsson
                errcode = 0 
            else
                h_id = 0
                errcode = c_io_error 
            end
            system.h_id = h_id
        elseif (lname=="t_end")
            system.t_end = parse(Float64, var) * units_fact
            errcode = 0
        elseif (lname=="p_total" || lname=="total_pressure")
            system.total_pressure = parse(Float64, var) * units_fact
            errcode = 0
        elseif (lname=="vsheath_solve_method")
            lvar = lowercase(var)
            if (lvar == "ohmic_power")
                solve_method = s_ohmic_power 
                errcode = 0
            elseif (lvar == "flux_balance")
                solve_method = s_flux_balance
                errcode = 0 
            elseif (lvar == "flux_interpolation")
                solve_method = s_flux_interpolation 
                errcode = 0 
            else
                solve_method = 0
                errcode = c_io_error 
            end
            system.Vsheath_solving_method = solve_method 
        elseif (lname=="prerun")
            system.prerun = parse(Bool, var)
            errcode = 0
        elseif (lname=="plasma_potential_guess")
            system.plasma_potential = parse(Float64, var)
            errcode = 0
        elseif (lname=="t_e_min")
            system.T_e_min = parse(Float64, var) * units_fact
            errcode = 0
        elseif (lname=="t_e_max")
            system.T_e_max = parse(Float64, var) * units_fact
            errcode = 0
        end
    else
        errcode = 0
    end
    return errcode 
end


function EndSystemBlock!(read_step::Int64, system::System)

    errcode = 0

    if (read_step == 1)
        if (system.A == 0)
            if (system.radius > 0 && system.l > 0)
                system.A = 2.0*pi*system.radius^2 +
                    2.0*pi*system.radius*system.l
            else
                print("***ERROR*** System area has not been defined\n")
                return c_io_error 
            end
        end
        if (system.V == 0)
            if (system.radius > 0 && system.l > 0)
                system.V = pi*system.radius^2 * system.l
            else
                print("***ERROR*** System volume has not been defined\n")
                return c_io_error 
            end
        end
        if (system.l == 0)
            print("***ERROR*** System length has not been defined\n")
            return c_io_error 
        end
        if (system.radius == 0)
            if (system.V == 0 || system.A == 0)
                print("***ERROR*** System radius has not been defined\n")
                return c_io_error 
            end
        end
        if (system.h_id == 0)
            print("***ERROR*** System h factor has not been defined\n")
            return c_io_error 
        end
        if (system.drivf == 0 && system.P_shape != p_constant )
            print("***ERROR*** System driving frequency has not been defined\n")
            return c_io_error 
        end
        if (system.drivP == 0)
            print("***ERROR*** System input power has not been defined\n")
        elseif system.drivP < 0
            print("***ERROR*** System input power must be positive\n")
        end

        if (system.P_duty_ratio < 0.0 || system.P_duty_ratio > 1.0)
            print("***ERROR*** Power duty ratio is out of range\n")
            return c_io_error 
        end

        if (system.dt_start < 0.0)
            print("***ERROR*** Initial simulation time-step must be positive\n")
            return c_io_error 
        end

        if (system.t_end == 0)
            print("***ERROR*** Simulation time must be > 0 \n")
            return c_io_error 
        end

        if (system.total_pressure > 0)
            print("***NOTE*** TOTAL PRESSURE has been declared in SYSTEM block\n")
            print("  - Species block must be declared AFTER the system block\n")
            print("  - Species pressures declarations (if declared) must be between 0 and 1\n")
            print("  - Output blocks with partial pressure parameters must be declared AFTER the system block\n")
        end

        if (system.T_e_min <= 0)
            print("***ERROR*** Make sure that T_e_min > 0\n")
            return c_io_error 
        end

        if (system.on_slope < 0.0) || (system.on_slope >= 1.0)
            print("***ERROR*** Power on-slope must be >= 0.0 and < 1.0\n")
            return c_io_error 
        end

        if (system.off_slope < 0.0) || (system.off_slope >= 1.0)
            print("***ERROR*** Power off-slope must be >= 0.0 and < 1.0\n")
            return c_io_error 
        end

        if (system.off_slope + system.on_slope >= 1.0)
            print("***ERROR*** Power on and off-slope cannot overlap\n")
            return c_io_error 
        end


        system.Lambda = GetLambda(system)
        system.drivOmega = system.drivf * 2.0 * pi
        system.P_absorbed = system.drivP / system.V
    end
    return errcode
end
    

function EndFile_System!(read_step::Int64, system::System)
    errcode = 0
    return errcode
end


function GetUnits!(var::Union{String,SubString{String}})
    # Identify units definition
    units_fact = 1.0
    units_index = findlast("_", var)
    if !(units_index === nothing)
        units_index = units_index[1]
        units_str = lowercase(var[units_index+1:end])
        match_flag = false
        if (units_str=="ev")
            # this eV refers to TEMPERATURE!
            units_fact = 1.0/K_to_eV
            match_flag = true
        elseif (units_str=="eev")
            # this eV refers to ENERGY 
            units_fact = e
            match_flag = true
        elseif (units_str=="mtorr")
            units_fact = 0.13332237
            match_flag = true
        elseif (units_str=="kw" || units_str=="khz")
            units_fact = 1.e3
            match_flag = true
        elseif (units_str=="mhz" || units_str=="mw")
            units_fact = 1.e6
            match_flag = true
        elseif (units_str=="ghz")
            units_fact = 1.e9
            match_flag = true
        elseif (units_str=="microns")
            units_fact = 1.e-6
            match_flag = true
        elseif (units_str=="nm")
            units_fact = 1.e-9
            match_flag = true
        elseif (units_str=="sccm")
            # The units_fact still needs to be divided by system.V, however,
            # just in case the volume changes, this is done in FunctionTerms.jl
            ns = 2.686780111798444e25 # Standard density at Ps = 101325 Pa and Ts = 273.15 K
            units_fact = 1.e-6 / 60 * ns
            match_flag = true
        end
        if match_flag
            var = strip(var[1:units_index-1])
        end
    end
    return units_fact, var
end

end