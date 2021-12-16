module InputBlock_System

using SharedData: c_io_error
using SharedData: p_ccp_id, p_icp_id
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation
using SharedData: System
using SharedData: K_to_eV
using PlasmaParameters: GetLambda


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

function StartFile_System!(read_step::Int64, system::System) 

    errcode = 0

    return errcode
end


function StartSystemBlock!(read_step::Int64, system::System)

    errcode = 0
    if (read_step == 1)
        system.A = 0.0
        system.V = 0.0
        system.l = 0.0
        system.radius = 0.0
        system.power_input_method = 0
        system.Vsheath_solving_method = 0
        system.drivf = 0.0
        system.drivOmega = 0.0
        system.drivP = 0.0
        system.t_end = 0.0
        system.Lambda = 0.0
        system.prerun = true
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
        elseif (lname=="electrode_area")
            system.electrode_area = parse(Float64,var)
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
        elseif (lname=="voltage")
            system.drivV = parse(Float64, var)
            errcode = 0
        elseif (lname=="power_method" || lname=="input_power_method" ||
                lname=="power_input_method")
            lvar = lowercase(var)
            if (lvar == "ccp")
                p_id = p_ccp_id
                errcode = 0
            elseif (lvar == "icp")
                p_id = p_icp_id
                errcode = 0 
            else
                p_id = 0
                errcode = c_io_error 
            end
            system.power_input_method = p_id
        elseif (name == "I" || lname == "current" ||
            lname == "driving_current")
            system.drivI = parse(Float64,var)
            errcode = 0
        elseif (lname=="t_end")
            system.t_end = parse(Float64, var) * units_fact
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
        if (system.power_input_method == 0)
            print("***ERROR*** System input power method has not been defined\n")
            return c_io_error 
        end
        if (system.drivf == 0)
            print("***ERROR*** System driving frequency has not been defined\n")
            return c_io_error 
        end
        if (system.drivP == 0)
            if (system.drivV == 0 && system.drivI == 0)
                print("***ERROR*** System input power has not been defined\n")
                return c_io_error 
            else
                system.drivP =  system.drivI * system.drivV
            end
        end

        if (system.t_end == 0)
            print("***ERROR*** Simulation time must be > 0 \n")
            return c_io_error 
        end

        system.Lambda = GetLambda(system)
        system.drivOmega = system.drivf * 2.0 * pi
    end
    return errcode
end
    

function EndFile_System!(read_step::Int64, system::System)
    errcode = 0
    return errcode
end


function GetUnits!(var::SubString{String})
    # Identify units definition
    units_fact = 1.0
    units_index = findlast("_", var)
    if !(units_index === nothing)
        units_index = units_index[1]
        units_str = lowercase(var[units_index+1:end])
        match_flag = false
        if (units_str=="ev")
            units_fact = 1.0/K_to_eV
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