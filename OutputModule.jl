module OutputModule

using SharedData: Species, Reaction, System, SpeciesID, Output
using SharedData: kb, K_to_eV
using PlasmaParameters: UpdateSpeciesParameters!
using PlasmaSheath: GetSheathVoltage
using WallFlux: UpdatePositiveFlux!, UpdateNegativeFlux!
using SolveSystem: ExecuteProblem

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
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID, output_list::Vector{Output})
    
    # Initial conditions
    dens0, temp0 = SetupInitialConditions(species_list)

    output_tuple_list = Tuple[]
    for output in output_list
        for output_flag in output.output_flag_list
            init = cat(temp0, dens0,dims=1)
            tspan = (0, system.t_end)
            GM_tuple = (system, sID, species_list, reaction_list)

            # The output flag come from a output-structure
            if output_flag == o_pL
                output_tuple = Output_PL(GM_tuple, temp0, dens0)

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
    dens::Vector{Float64})

    sID = GM_tuple[2]
    system = GM_tuple[1]

    TAr = temp[sID.Ar]
    Te = Float64[]
    ne = Float64[]
    pL = Float64[]
    n_steps = 150
    for e in [10^y for y in range(log10(19.5), log10(25), length=n_steps)]
        dens[sID.Ar] = 9.5 * 10^e
        init = cat(temp,dens,dims=1)
        p = dens[sID.Ar] * kb * TAr
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
        push!(Te, sol[1,end]*K_to_eV)
        push!(ne, sol[4,end])
    end

    return ("pL", Te, ne, pL)
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


end