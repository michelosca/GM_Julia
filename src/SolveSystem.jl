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

module SolveSystem

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: e, K_to_eV
using SharedData: o_single_run, c_io_error
using SharedData: p_constant, p_square
using PlasmaParameters: UpdateParameters!
using PlasmaSheath: GetSheathVoltage!
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using FunctionTerms: GetDensRateFunction, GetTempRateFunction
using DifferentialEquations: ODEProblem, solve, Trapezoid, Rosenbrock23, ImplicitEuler
using DifferentialEquations: DiscreteCallback, VectorContinuousCallback
using DifferentialEquations: CallbackSet
using DifferentialEquations: terminate!, set_proposed_dt!, get_proposed_dt 
using Printf
using PrintModule: PrintErrorMessage, PrintWarningMessage, PrintMessage
#using PrintModule: PrintSimulationState

function ExecuteProblem(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID::SpeciesID,
    o_case::Int64=0)

    # Initial conditions
    dens, temp = GetInitialConditions(species_list)
    init = cat(temp,dens,dims=1)
    tspan = (0, system.t_end)
    p = (system, sID, species_list, reaction_list )

    # Event handling
    cb = nothing
    if system.P_shape == p_square
        cb_duty_ratio = VectorContinuousCallback(condition_pulsed_power,
            affect_pulsed_power!, 2, affect_neg! = nothing, save_positions=(true,true))
        cb_error = DiscreteCallback(condition_error, affect_error!)
        cb_Te_lower_bound = DiscreteCallback(condition_Te_lower_bound, affect_Te_lower_bound!)
        cb = CallbackSet(cb_duty_ratio, cb_Te_lower_bound, cb_error)
    else
        cb_error = DiscreteCallback(condition_error, affect_error!)
        cb = CallbackSet(cb_error)
    end

    # ODE problem
    prob = ODEProblem{true}(ode_fn!, init, tspan, p)
    if o_case == o_single_run
        save_flag = true
    else
        save_flag = false
    end

    PrintMessage(system, "Running problem ...\n")
    #PrintSimulationState(temp, dens, species_list, system, sID)
    sol = solve(prob,
        Rosenbrock23(autodiff=false),
        dt = system.dt_start,
        #abstol = 1.e-8,
        #reltol = 1.e-6,
        maxiters = 1.e7,
        callback = cb,
        save_everystep = save_flag
    )
    PrintMessage(system, "Problem execution finished\n")
    return sol
end

function ode_fn!(dy::Vector{Float64}, y::Vector{Float64}, p::Tuple, t::Float64)
    # dy are the derivatives of y, i.e. the lhs of the ODE equations
    # y are [temperature,dens] (the input Vector in other functions) 
    # p are the free parameters the y-functions may have
    # t is time

    # Make sure lenght of y is the same as p
    system = p[1]
    sID = p[2]
    species_list = p[3]
    reaction_list = p[4]
    n_species = length(species_list)
    dens = y[2:end]
    temp = Float64[] 
    for s in species_list
        if s.id== sID.electron
            push!(temp, y[1])
        else
            push!(temp, s.temp)
        end
    end

    # Update species parameters
    errcode = UpdateParameters!(temp, dens, species_list, reaction_list, system, sID)
    if errcode == c_io_error
        system.errcode = errcode 
        PrintErrorMessage(system, "UpdateSpeciesParameters failed")
    end

    # First get the positive ion fluxes
    errcode = UpdatePositiveFlux!(species_list)
    if errcode == c_io_error
        system.errcode = errcode 
        PrintErrorMessage(system, "UpdatePositiveFlux failed")
    end

    # Calculate the sheath potential
    errcode = GetSheathVoltage!(system, species_list, sID, t)
    if errcode == c_io_error
        system.errcode = errcode 
        PrintErrorMessage(system, "GetSheathVoltage failed")
    end

    # Calculate the electron flux  
    errcode = UpdateNegativeFlux!(species_list, system, sID)
    if errcode == c_io_error
        system.errcode = errcode 
        PrintErrorMessage(system, "UpdateNegativeFlux failed")
    end

    ### Generate the dy array
    # Temperature equation
    dy[1] = GetTempRateFunction(temp, dens, species_list[sID.electron],
        species_list, reaction_list, system, sID, t)
    for i in 1:n_species
        #id = i - ((i-1)÷n_species)*n_species
        species = species_list[i]

        # Density equation
        dy[i+1] = GetDensRateFunction(temp, dens, species,
            species_list, reaction_list, system, sID)
    end
end


function GetInitialConditions(species_list::Vector{Species})
    dens = Float64[]
    temp = Float64[]
    for s in species_list
        push!(dens, s.dens)
        if s.has_temp_eq
            push!(temp, s.temp)
        end
    end
    return dens, temp
end


function condition_pulsed_power(out, u, t, integrator)
    # Event when time has past duty cycle 
    p = integrator.p
    system = p[1]

    t_start = system.P_start
    t_offset = t - t_start
    freq = system.drivf
    dr = system.P_duty_ratio
    if t_offset > 0
        period_fraction = t_offset * freq
        out[1] = period_fraction - floor(period_fraction) - dr
        out[2] = period_fraction - round(period_fraction)
    else
        out[1] = -dr
        out[2] = 0.0
    end

    #T_e_min = system.T_e_min
    #T_e = integrator.u[1]
    #out[3] = T_e_min - T_e
end


function affect_pulsed_power!(integrator, cb_index)
    # What to do when the event occurs
    p = integrator.p
    system = p[1]
    dt = system.dt_start
    if cb_index == 1
        # Power off
        system.P_absorbed = 0.0 
        #message = @sprintf(" Power switch off. Time = %10g s; dt = %10g s\n", integrator.t, dt)
        #PrintMessage(system, message)
    elseif cb_index == 2
        # Power on
        system.P_absorbed = system.drivP / system.V 
        #message = @sprintf(" Power switch on.  Time = %10g s; dt = %10g s\n", integrator.t, dt)
        #PrintMessage(system, message)
    #elseif cb_index == 3
    #    T_e = integrator.u[1]
    #    T_e_min = system.T_e_min
    #    if T_e < T_e_min
    #        integrator.u[1] = system.T_e_min
    #    end
    end
    set_proposed_dt!(integrator, dt) 
end


function condition_error(u, t, integrator)
    # Event when time has past duty cycle 
    p = integrator.p
    system = p[1]
    return system.errcode == c_io_error
end


function affect_error!(integrator)
    p = integrator.p
    system = p[1]
    PrintErrorMessage(system, "Abort simulation")
    terminate!(integrator)
end


function condition_Te_lower_bound(u, t, integrator)
    u = integrator.u
    system = integrator.p[1]

    if u[1] < system.T_e_min
        return true 
    end

    return false 
end


function affect_Te_lower_bound!(integrator)
    p = integrator.p
    system = p[1]
    integrator.u[1] = system.T_e_min
end

end
