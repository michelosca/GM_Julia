module OutputModule

using SharedData: Species, Reaction, System, Output
using SolveSystem: ExecuteProblem
using SharedData: kb, K_to_eV
using InputBlock_Species: s_electron_id, s_Ar_id
using PlasmaSheath: GetSheathVoltage
using WallFlux: GetIonFlux, GetElectronFlux!
using GenerateODEs: temp_eq_elastic, temp_eq_gainloss
using GenerateODEs: temp_eq_ethreshold, temp_eq_flux, temp_eq_inpower

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
const o_pL = 1
const o_singel_sol = 2

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GenerateOutputs
#   - SetupInitialConditions
#   - Output_PL

function GenerateOutputs(
    dens_eq::Vector{Any}, temp_eq::Vector{Any},
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system_list::Vector{System}, output_list::Vector{Output})
    
    # Join dens and temp equation list into one list
    eq_list = cat(temp_eq, dens_eq,dims=1)

    # Initial conditions
    dens0, temp0 = SetupInitialConditions(species_list)

    output_tuple_list = Tuple[]
    for system in system_list
        for output in output_list
            for output_flag in output.output_flag_list
                init = cat(temp0, dens0,dims=1)
                tspan = (0, system.t_end)
                GM_tuple = (eq_list, system, species_list, reaction_list)

                # The output flag come from a output-structure
                if output_flag == o_pL
                    output_tuple = Output_PL(GM_tuple, temp0, dens0, system,
                        temp_eq, species_list, reaction_list)

                elseif (output_flag == o_singel_sol)
                    print("Solving single problem...\n")
                    sol = ExecuteProblem(init, tspan, GM_tuple)
                    output_tuple = ("single sol", sol)
                end

                push!(output_tuple_list, output_tuple)

            end # Output flags loop
        end # Output list loop
    end # System list loop

    return output_tuple_list
end

function Output_PL(GM_tuple::Tuple, temp::Vector{Float64},
    dens::Vector{Float64}, system::System, temp_eq::Vector{Any}, species_list,
    reaction_list)

    TAr = temp[s_Ar_id]
    Te = Float64[]
    ne = Float64[]
    pL = Float64[]
    n_terms = length(temp_eq[s_electron_id])
    n_steps = 200
    energy_terms = zeros(n_steps, n_terms)
    istep = 1
    for e in [10^y for y in range(log10(19.5), log10(24), length=n_steps)]
        dens[s_Ar_id] = 9.5 * 10^e
        init = cat(temp,dens,dims=1)
        p = dens[s_Ar_id] * kb * TAr
        push!(pL, p * system.l)
        print("P*L = ",p*system.l,";Executing GM problem...\n")
        if (pL[end] > 1.e4) # n = 1.e25
            t_end = 1.0 
        elseif (pL[end] > 1.e3) # n = 1.e24
            t_end = 3.e-1 
        elseif (pL[end] > 1.e2) # n = 1.e23
            t_end = 2.e-1 # pL is between 1.e3 and 1.e2
        elseif (pL[end] > 1.e1) # n = 1.e22
            t_end = 1.e-1 # pL is between 1.e2 and 1.e1
        elseif (pL[end] > 1.e0) # n = 1.e21
            t_end = 2.e-2 # pL is between 1.e1 and 1.e0
        else
            t_end = 2.e-3 # pL is lower than 1.e0
        end
        tspan = (0,t_end)
        sol = ExecuteProblem(init, tspan, GM_tuple)
        term_list = GetTempEquationTerms(temp_eq[s_electron_id], dens, temp, species_list,
            reaction_list, system)
        push!(Te, sol[1,end]*K_to_eV)
        push!(ne, sol[5,end])
        energy_terms[istep,:] = term_list
        istep += 1
    end

    return ("pL", Te, ne, energy_terms, pL)
end


function SetupInitialConditions(species_list::Vector{Species})
    dens = Float64[]
    temp = Float64[]
    for s in species_list
        push!(dens, s.dens0)
        push!(temp, s.temp0)
    end
    return dens, temp
end

function GetTempEquationTerms(temp_eq::Vector{Tuple}, dens, temp, species_list,
    reaction_list, system)

    # First get the positive ion fluxes
    flux_list = GetIonFlux(dens, temp, species_list, reaction_list, system)

    # Calculate the sheath potential
    V_sheath = GetSheathVoltage(dens, temp, species_list, reaction_list,
        system, flux_list)

    # Calculate the electron flux  
    GetElectronFlux!(flux_list, species_list)

    value_list = Float64[]
    for current_tuple in temp_eq 
        flag = current_tuple[1]
        funct = current_tuple[2]

        # Temperature equations
        if (flag == temp_eq_gainloss)
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
            f_curr = funct(dens, temp, species_list, flux_list, system, V_sheath)
        elseif (flag == temp_eq_inpower)
            # Power absorption
            f_curr = funct(dens, system)
        else
            f_curr = 0.0
            print("***WARNING*** Temperature equation term not found!\n")
        end
        push!(value_list, f_curr)
    end
    return value_list
end

end