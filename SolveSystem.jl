module SolveSystem

using SharedData: System, Species
using SharedData: e#, K_to_eV
using PlasmaParameters: UpdateSpeciesParameters!
using PlasmaSheath: GetSheathVoltage
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using FunctionTerms: GetDensRateFunction, GetTempRateFunction
using DifferentialEquations: ODEProblem, solve, Tsit5, Euler 

function ExecuteProblem(init::Vector{Float64}, tspan::Tuple, p::Tuple)

    prob = ODEProblem(ode_fn!, init, tspan, p)
    sol = solve(prob, Tsit5(), maxiters=1.e7, reltol=1e-8, abstol=1e-8)#, dt = 1.e-10)
    return sol
end

function ode_fn!(dy::Vector{Float64}, y::Vector{Float64}, p::Tuple, t::Float64)
    # dy are the derivatives of y, i.e. the lhs of the ODE equations
    # y are [temperature,dens] (the input Vector in other functions) 
    # p are the free parameters the y-functions may have
    # t is time

    #print("Time ", t,"\n")

    # Make sure lenght of y is the same as p
    system = p[1]
    sID = p[2]
    species_list = p[3]
    reaction_list = p[4]
    n_species = length(species_list)
    dens = y[n_species + 1:end]
    temp = y[1:n_species]
    #print(t, " densities ", dens[[1,3,4,6]],"\n neutrality ", dens[1]+dens[4]-dens[3]-dens[6],"\n")

    #Te_eV = temp[sID.electron] * K_to_eV
    #if Te_eV > 7 
    #    print("***WARNING*** Te above threshold ", Te_eV," eV @ t =", t,"\n")
    #end
    # Update species parameters
    UpdateSpeciesParameters!(temp, dens, species_list, system, sID)

    # First get the positive ion fluxes
    UpdatePositiveFlux!(species_list, system)

    # Calculate the sheath potential
    V_sheath = GetSheathVoltage(species_list, system, sID)

    # Calculate the electron flux  
    UpdateNegativeFlux!(species_list, V_sheath)
    
    # Generate the dy array
    #print("Start\n")
    #print(" - Time ", t,"\n")
    #print(" - V_sheat ", V_sheath,"\n")
    dens_balance = 0.0
    flux_balance = 0.0
    for i in 1:n_species
        id = i - ((i-1)÷n_species)*n_species
        species = species_list[id]

        # Temperature equation
        dy[i] = GetTempRateFunction(temp, dens, species, species_list,
            reaction_list, system, V_sheath, sID)
        # Density equation
        dy[i+n_species] = GetDensRateFunction(temp, dens, species,
            reaction_list, system, sID)

        if (!(species.charge==0))
            dens_balance += species.dens * species.charge
            flux_balance += species.flux * species.charge
            #print(" - Species ", species.id, "\n")
            #print("   - Species dens: ", species.dens, " temp ", species.temp,"\n")
            #print("   - Species flux: ", species.flux,"\n")
            #print("   - Temp FUNCTION: ",dy[i],"\n")
            #print("   - Dens FUNCTION: ",dy[i+n_species],"\n")
        end
    end
    dens_balance /= e
    flux_balance /= e
    #print(t," - Dens BALANCE: ",dens_balance)
    #print(" - Flux BALANCE: ",flux_balance)
    #print(" - V_sheath: ",V_sheath,"\n")
    #print("\n\n")
end

end