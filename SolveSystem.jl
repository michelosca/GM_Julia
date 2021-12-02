module SolveSystem

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: e, K_to_eV
using PlasmaParameters: UpdateSpeciesParameters!
using PlasmaSheath: GetSheathVoltage
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using FunctionTerms: GetDensRateFunction, GetTempRateFunction
using DifferentialEquations: ODEProblem, solve, Tsit5, Euler 

function ExecuteProblem(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID::SpeciesID)

    # Initial conditions
    dens, temp = GetInitialConditions(species_list)
    init = cat(temp,dens,dims=1)

    tspan = (0, system.t_end)
    p = (system, sID, species_list, reaction_list )


    prob = ODEProblem(ode_fn!, init, tspan, p)

    print("Solving single problem...\n")
    sol = solve(prob, Tsit5(), maxiters=1.e7, dt = 1.e-14)#, dt = 1.e-10, maxiters=1.e7) #, reltol=1e-8, abstol=1e-8, dt = 1.e-13)

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
    dens = y[n_species + 1:end]
    temp = y[1:n_species]

    # Update species parameters
    UpdateSpeciesParameters!(temp, dens, species_list, system, sID)

    # First get the positive ion fluxes
    UpdatePositiveFlux!(species_list, system)

    # Calculate the sheath potential
    V_sheath = GetSheathVoltage(species_list, system, sID)

    # Calculate the electron flux  
    UpdateNegativeFlux!(species_list, system, sID, V_sheath)
    
    # Generate the dy array
    for i in 1:n_species
        id = i - ((i-1)÷n_species)*n_species
        species = species_list[id]

        # Temperature equation
        dy[i] = GetTempRateFunction(temp, dens, species, species_list,
            reaction_list, system, V_sheath, sID)
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