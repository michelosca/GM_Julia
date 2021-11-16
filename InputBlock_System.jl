module InputBlock_System

using SharedData: c_io_error
using SharedData: p_ccp_id, p_icp_id
using SharedData: System

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
global system_list = System[]

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - StartFile_System
# - StartSystemBlock
# - ReadSystemEntry
# - EndSystemBlock
#   - AddSystemToList

function StartFile_System(read_step) 

    errcode = 0

    if (read_step == 1)
        global system_list = System[]
    end

    return errcode
end


function StartSystemBlock(read_step::Int64)

    errcode = 0

    if (read_step == 1)
        global input_area = 0.0
        global input_vol = 0.0
        global input_len = 0.0
        global input_radius = 0.0
        
        global input_power_method = 0 
        global input_drivf = 0.0
        global input_drivP = 0.0
        global input_drivV = 0.0
        global input_drivI = 0.0
        global input_time = 0.0

        global input_electrodeA = 0.0 
    end

    return errcode
end


function ReadSystemEntry(name, var, read_step)

    errcode = c_io_error 

    if (read_step == 1)
        # Identify units definition
        units_fact = 1.0
        units_index = findlast("_", var)
        if (units_index === nothing)
            units = ""
        else
            units_index = units_index[1]
            units_str = lowercase(var[units_index+1:end])
            if (units_str=="kw" || units_str=="khz")
                units_fact = 1.e3
                var = var[1:units_index-1]
            elseif (units_str=="mhz" || units_str=="mw")
                units_fact = 1.e6
                var = var[1:units_index-1]
            elseif (units_str=="ghz")
                units_fact = 1.e9
                var = var[1:units_index-1]
            elseif (units_str=="microns")
                units_fact = 1.e-6
                var = var[1:units_index-1]
            end
        end
        
        # Identify variable
        lname = lowercase(name)
        if (name=="A" || lname=="area")
            global input_area = parse(Float64,var)
            errcode = 0
        elseif (lname=="electrode_area")
            global input_electrodeA = parse(Float64,var)
            errcode = 0
        elseif (lname=="v" || lname=="volume" || lname=="vol")
            global input_vol = parse(Float64, var)
            errcode = 0
        elseif (lname=="l" || lname=="length")
            global input_len = parse(Float64, var) 
            errcode = 0
        elseif (lname=="radius" || lname=="r")
            global input_radius = parse(Float64, var) 
            errcode = 0
        elseif (name=="f" || lname=="frequency" ||
            lname=="freq" || lname=="driving_frequency")
            global input_drivf = parse(Float64,var) * units_fact
            errcode = 0
        elseif (name=="P" || lname=="power_input" ||
            lname=="power")
            global input_drivP = parse(Float64, var) * units_fact
            errcode = 0
        elseif (lname=="voltage")
            global input_drivV = parse(Float64, var)
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
            global input_power_method = p_id
        elseif (name == "I" || lname == "current" ||
            lname == "driving_current")
            global drivI = parse(Float64,var)
            errcode = 0
        elseif (lname=="t_end")
            global input_time = parse(Float64, var) * units_fact
            errcode = 0
        end
    else
        errcode = 0
    end
    return errcode 
end


function EndSystemBlock(read_step::Int64)

    errcode = 0

    if (read_step == 1)
        if (input_area == 0)
            if (input_radius > 0)
                global input_area = pi * input_radius^2
            else
                print("***ERROR*** System area has not been defined\n")
                return c_io_error 
            end
        end
        if (input_vol == 0)
            if (input_radius > 0 && input_len > 0)
                global input_vol = 2.0 * pi * input_radius * input_len
            else
                print("***ERROR*** System volume has not been defined\n")
                return c_io_error 
            end
        end
        if (input_len == 0)
            print("***ERROR*** System length has not been defined\n")
            return c_io_error 
        end
        if (input_radius == 0)
            if (input_vol == 0 || input_area == 0)
                print("***ERROR*** System radius has not been defined\n")
                return c_io_error 
            end
        end
        if (input_power_method == 0)
            print("***ERROR*** System input power method has not been defined\n")
            return c_io_error 
        elseif (input_power_method == p_ccp_id)
            if (input_electrodeA == 0.0)
                print("***ERROR*** The 'electrode_area' has not been defined\n")
            end
        end
        if (input_drivf == 0)
            print("***ERROR*** System driving frequency has not been defined\n")
            return c_io_error 
        end
        if (input_drivP == 0)
            if (input_drivV == 0 || input_drivI == 0)
                print("***ERROR*** System input power has not been defined\n")
                return c_io_error 
            else
                global input_drivP =  input_drivI * input_drivV
            end
        end
        if (input_drivV == 0)
            if (input_drivP == 0 || input_drivI == 0)
                print("***ERROR*** System input voltage has not been defined\n")
                return c_io_error 
            else
                global input_drivV =  input_drivP / input_drivI
            end
        end
        if (input_drivI == 0)
            if (input_drivP == 0 || input_drivV == 0)
                print("***ERROR*** System input current has not been defined\n")
                return c_io_error 
            else
                global input_drivI =  input_drivP / input_drivV
            end
        end
        if (input_time == 0)
            print("***ERROR*** Simulation time must be > 0 \n")
            return c_io_error 
        end

        errcode = AddSystemToList(input_area, input_vol, input_len,
            input_power_method, input_drivf, input_drivP, input_drivV,
            input_drivI, input_time, input_radius, input_electrodeA)
    end
    return errcode
end

function AddSystemToList(area::Float64, vol::Float64, l::Float64,
    p_method::Int64, drivf::Float64, drivP::Float64,
    drivV::Float64, drivI::Float64, time::Float64, radius::Float64,
    electrodeA::Float64)
    
    drivOmega = drivf * 2.0 * pi

    errcode = 0
    try
        # Add react to reaction_list
        system = System(area, vol, l, radius, p_method, drivf, drivOmega, drivP,
            drivV, drivI, time, electrodeA)
        push!(system_list, system)
    catch
        print("***ERROR*** While attaching system\n")
        errcode = c_io_error 
    end
    return errcode
end

end