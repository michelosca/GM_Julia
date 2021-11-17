module OutputModule

using SharedData: Species, Reaction, System, Output
using InputBlock_Species: s_electron_id
using SolveSystem: ExecuteProblem
using SharedData: kb, K_to_eV
using InputBlock_Species: s_electron_id, s_Ar_id
using PlasmaSheath: GetSheathVoltage
using PlasmaParameters: UpdateSpeciesParameters!
using WallFlux: UpdateIonFlux!, UpdateElectronFlux!
using GenerateODEs: temp_eq_elastic, temp_eq_gainloss
using GenerateODEs: temp_eq_ethreshold, temp_eq_flux, temp_eq_inpower
using GenerateODEs: dens_eq_gainloss, dens_eq_flux, eq_empty

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
    system::System, output_list::Vector{Output})
    
    # Join dens and temp equation list into one list
    eq_list = cat(temp_eq, dens_eq,dims=1)

    # Initial conditions
    dens0, temp0 = SetupInitialConditions(species_list)

    output_tuple_list = Tuple[]
    for output in output_list
        for output_flag in output.output_flag_list
            init = cat(temp0, dens0,dims=1)
            tspan = (0, system.t_end)
            GM_tuple = (eq_list, system, species_list, reaction_list)

            # The output flag come from a output-structure
            if output_flag == o_pL
                output_tuple = Output_PL(GM_tuple, temp0, dens0, system,
                    temp_eq, species_list)

            elseif (output_flag == o_singel_sol)
                print("Solving single problem...\n")
                sol = ExecuteProblem(init, tspan, GM_tuple)
                output_tuple = ("single sol", sol)
            end

            push!(output_tuple_list, output_tuple)

        end # Output flags loop
    end # Output list loop

    return output_tuple_list
end

function Output_PL(GM_tuple::Tuple, temp::Vector{Float64},
    dens::Vector{Float64}, system::System, temp_eq::Vector{Any},
    species_list::Vector{Species})

    TAr = temp[s_Ar_id]
    Te = Float64[]
    ne = Float64[]
    pL = Float64[]
    n_terms = length(temp_eq[s_electron_id])
    n_steps = 100
    energy_terms = zeros(n_steps, n_terms)
    istep = 1
    for e in [10^y for y in range(log10(19.5), log10(22), length=n_steps)]
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
        term_list = GetTempEquationTerms(dens, temp, s_electron_id, species_list,
            system, temp_eq[s_electron_id])
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
        push!(dens, s.dens)
        push!(temp, s.temp)
    end
    return dens, temp
end

function GetTempEquationTerms(dens::Vector{Float64}, temp::Vector{Float64},
    species_id::Int64, species_list::Vector{Species}, system::System,
    funct_list::Vector{Tuple})

    UpdateSpeciesParameters!(temp, dens, species_list, system)


    # First get the positive ion fluxes
    UpdateIonFlux!(species_list, system)

    # Calculate the sheath potential
    V_sheath = GetSheathVoltage(species_list, system)

    # Calculate the electron flux  
    UpdateElectronFlux!(species_list)

    value_list = Float64[]
    for current_tuple in funct_list
        flag = current_tuple[1]
        funct = current_tuple[2]

        # Density equations
        if (flag == dens_eq_gainloss)
            # Particle gain/loss
            f_curr = funct(dens,temp)
        elseif (flag == dens_eq_flux)
            # Particle fluxes
            f_curr = funct(species_list[species_id])

        # Temperature equations
        elseif (flag == temp_eq_gainloss || flag == temp_eq_elastic ||
            flag == temp_eq_ethreshold)
            # Particle gain/loss, elastic or ethreshold
            f_curr = funct(dens,temp)
        elseif (flag == temp_eq_flux)
            # Energy fluxes
            f_curr = funct(dens, temp, species_list, V_sheath)
        elseif (flag == temp_eq_inpower)
            # Power absorption
            f_curr = funct(dens, system)
        elseif (flag == eq_empty)
            continue
        end
        push!(value_list, f_curr)
    end

    return value_list
end

end