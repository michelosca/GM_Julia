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
using PlasmaParameters: UpdateSpeciesParameters!
using PlasmaSheath: GetSheathVoltage!
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using FunctionTerms: GetDensRateFunction, GetTempRateFunction
using DifferentialEquations: ODEProblem, solve, Trapezoid, Rodas5, Rosenbrock23
using Printf
using PrintModule: PrintErrorMessage, PrintWarningMessage

function ExecuteProblem(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID::SpeciesID,
    o_case::Int64=0)

    # Initial conditions
    dens, temp = GetInitialConditions(species_list)
    init = cat(temp,dens,dims=1)

    tspan = (0, system.t_end)
    p = (system, sID, species_list, reaction_list )


    prob = ODEProblem{true}(ode_fn!, init, tspan, p)

    if (o_case == o_single_run)
        save_flag = true
    else
        save_flag = false
    end

    try
        print("Solving single problem ...\n")
        #PrintSimulationState(temp, dens, species_list, system, sID)
        sol = solve(prob,
            Trapezoid(autodiff=false),
            dt=1.e-12,
            #abstol=1.e-10,
            #reltol=1.e-6,
            maxiters=1.e7,
            save_everystep=save_flag
        )
        return sol
    catch
        PrintWarningMessage(system, "Re-run problem with increases solving tolerances")
        sol = solve(prob,
            #Trapezoid(autodiff=false),
            Rosenbrock23(autodiff=false),
            dt=1.e-12,
            abstol=1.e-10,
            reltol=1.e-6,
            maxiters=1.e7,
            save_everystep=save_flag
        )
        return sol
    end

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
    errcode = UpdateSpeciesParameters!(temp, dens, species_list, reaction_list, system, sID)
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


end