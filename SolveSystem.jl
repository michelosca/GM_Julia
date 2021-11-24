module SolveSystem

using SharedData: System, Species
using SharedData: e, K_to_eV
using PlasmaParameters: UpdateSpeciesParameters!
using PlasmaSheath: GetSheathVoltage
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using FunctionTerms: GetDensRateFunction, GetTempRateFunction
using DifferentialEquations: ODEProblem, solve, Tsit5, Euler 

function ExecuteProblem(init::Vector{Float64}, tspan::Tuple, p::Tuple)

    sol_list = []
    rate = 1.0

    system = p[1]
    sID = p[2]
    species_list = p[3]
    reaction_list = p[4]
    n_species = length(species_list)
    while rate > 0.001

        prob = ODEProblem(ode_fn!, init, tspan, p)
        if rate >= 1
            sol = solve(prob, Tsit5(), maxiters=1.e7, reltol=1e-8, abstol=1e-8, dt = 1.e-13)
        else
            sol = solve(prob, Tsit5(), maxiters=1.e7, reltol=1e-8, abstol=1e-8, dt = 1.e-10)
        end

        push!(sol_list, sol)

        # Update species parameters
        dens_0 = init[n_species + sID.electron]
        for s in species_list
            s.dens = sol[n_species + s.id, end]
            s.temp = sol[s.id, end]
            init[s.id] = s.temp
            init[s.id + n_species] = s.dens
        end
        p = (system, sID, species_list, reaction_list)

        # Check convergence rate
        dens_end = species_list[sID.electron].dens
        rate = abs(dens_end - dens_0)/dens_0
        print("Convergence rate: ", rate,"\n")

    end
    return sol_list
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
            reaction_list, system, sID)
    end
end

end