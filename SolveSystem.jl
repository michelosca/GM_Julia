module SolveSystem

using DifferentialEquations: ODEProblem, solve
using DifferentialEquations: Tsit5 
using WallFlux: MeanSheathVoltage, InterpolateSheathVoltage
using WallFlux: GetParticleFlux
using SharedData: p_ccp_id, p_icp_id 
using SharedData: System
using GenerateODEs: dens_eq_gainloss, dens_eq_flux, eq_empty
using GenerateODEs: temp_eq_elastic, temp_eq_gainloss
using GenerateODEs: temp_eq_ethreshold, temp_eq_flux, temp_eq_inpower

function ExecuteProblem(init::Vector{Float64}, tspan::Tuple, p::Tuple)

    prob = ODEProblem(ode_fn!, init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    return sol
end

function ode_fn!(dy,y,p,t)
    # dy are the derivatives of y, i.e. the lhs of the ODE equations
    # y are [temperature,dens] (the input Vector in other functions) 
    # p are the free parameters the y-functions may have
    # t is time

    # Make sure lenght of y is the same as p
    eq_list = p[1]
    system = p[2]
    species_list = p[3]
    reaction_list = p[4]
    n_species = length(species_list)
    temp = y[1:n_species]
    dens = y[n_species + 1:end]

    # Calculate the sheath potential
    if (system.power_input_method == p_ccp_id)
        V_sheath = MeanSheathVoltage(dens, system)
    elseif (system.power_input_method == p_ccp_id)
        V_sheath = InterpolateSheathVoltage(dens, temp, species_list,
            reaction_list, system)
    end

    # Calculate the particle fluxes
    flux_list = GetParticleFlux(dens, temp, species_list,
        reaction_list, system, V_sheath)
    #print("Fluxes ",flux_list,"\n")

    # Generate the dy array
    n = length(eq_list)
    for i in 1:n
        dy[i] = GatherListOfFunctions(dens, temp, flux_list, system, V_sheath,
            eq_list[i])
    end
end


function GatherListOfFunctions(dens::Vector{Float64}, temp::Vector{Float64},
    flux_list::Vector{Float64}, system::System, V_sheath::Float64,
    funct_list::Vector{Tuple})

    f_out = 0
    for current_tuple in funct_list
        flag = current_tuple[1]
        funct = current_tuple[2]

        # Density equations
        if (flag == dens_eq_gainloss)
            # Particle gain/loss
            f_curr = funct(dens,temp)
        elseif (flag == dens_eq_flux)
            # Particle fluxes 
            f_curr = funct(flux_list, system)

        # Temperature equations
        elseif (flag == temp_eq_gainloss)
            # Particle gain/loss
            f_curr = funct(dens,temp)
        elseif (flag == temp_eq_elastic)
            # Elastic collisions
            f_curr = funct(dens,temp)
        elseif (flag == temp_eq_ethreshold)
            # Intrinsic collision energy loss
            f_curr = funct(dens,temp)
        elseif (flag == temp_eq_flux)
            # Energy fluxes
            f_curr = funct(dens,temp, flux_list, system, V_sheath)
        elseif (flag == temp_eq_inpower)
            # Power absorption
            f_curr = funct(dens, system)
        elseif (flag == eq_empty)
            continue
        else
            print("Function flag is not recognized\n")
        end
        f_out += f_curr
    end

    return f_out
end

end