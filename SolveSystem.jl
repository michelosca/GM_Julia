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
using SharedData: o_single_run
using PlasmaParameters: UpdateSpeciesParameters!
using PlasmaSheath: GetSheathVoltage
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using FunctionTerms: GetDensRateFunction, GetTempRateFunction
using DifferentialEquations: ODEProblem, solve, Trapezoid, Rodas5
using Printf

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
        sol = solve(prob,
            Trapezoid(autodiff=false),
            dt=1.e-12,
            #abstol=1.e-8,
            #reltol=1.e-3,
            save_everystep=save_flag)
        return sol
    catch
        print("***WARNING*** Trapezoid solver failed. Rerunning with Rodas5\n")
        sol = solve(prob, Rodas5(autodiff=false), dt=1.e-12, abstol=1.e-8,
            reltol=1.e-3, save_everystep=save_flag)
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
    dens = y[n_species + 1:end]
    temp = y[1:n_species]

    # Update species parameters
    UpdateSpeciesParameters!(temp, dens, species_list, system, sID)

    # First get the positive ion fluxes
    UpdatePositiveFlux!(species_list, system)

    # Calculate the sheath potential
    GetSheathVoltage(species_list, system, sID, t)

    # Calculate the electron flux  
    UpdateNegativeFlux!(species_list, system, sID)

    # Generate the dy array
    for i in 1:n_species
        id = i - ((i-1)÷n_species)*n_species
        species = species_list[id]

        # Temperature equation
        dy[i] = GetTempRateFunction(temp, dens, species, species_list,
            reaction_list, system, sID, t)
        # Density equation
        dy[i+n_species] = GetDensRateFunction(temp, dens, species,
            species_list, reaction_list, system, sID)
    end
end


function GetInitialConditions(species_list::Vector{Species})
    dens = Float64[]
    temp = Float64[]
    for s in species_list
        push!(dens, s.dens)
        push!(temp, s.temp)
    end
    return dens, temp
end


end