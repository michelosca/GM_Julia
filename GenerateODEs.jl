module GenerateODEs

using SharedData: System, Species, Reaction
using SharedData: kb, e, eps0
using SharedData: p_ccp_id 
using InputBlock_Species: s_electron_id
using InputBlock_Reactions: r_elastic_id, r_excitat_id
using PowerInput: PowerInputFunction

###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - GenerateDensRateFunctionList
#   - GenerateDensRateFunction
#     - DensWallFluxFunction
#       - WallFluxFunction_CCP
#         - WallFluxFunction_CCP_ion_species
# - GenerateTempRateFunctionList
#   - GenerateTempRateFunction
#     - TempWallFluxFunction
#       - MeanSheathVoltage
#       - WallFluxFunction_CCP
#         - WallFluxFunction_CCP_ion_species
#     - PowerInputFunction


function GenerateDensRateFunctionList(system::System,
    species_list::Vector{Species}, reaction_list::Vector{Reaction})
    # Gathers the Density rate functions from GenerateDensRateFunction for 
    #each species in a list
    
    # Loop over all species
    dens_rate_funct_list = []
    for s in species_list
        sdens_funct_list = GenerateDensRateFunction(s, species_list,
            reaction_list, system)
        if (length(sdens_funct_list)>0)
            push!(dens_rate_funct_list, sdens_funct_list)
        end
    end
    return dens_rate_funct_list
end


function GenerateDensRateFunction(s::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction},
    system::System)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    sdens_funct_list = Function[]

    # Check if species is meant to have a density ODE 
    if (s.has_dens_eq)
        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s 
            if (r.id == r_excitat_id)
                continue
            end
            s_index = findall( x -> x == s.id, r.involved_species )
            if !(length(s_index) == 0)
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                # Only a density gain/loss term is added if there is a non-zero
                # particle balance
                if !(sign == 0)
                    push!(sdens_funct_list,
                        function (dens::Vector{Float64}, temp::Vector{Float64})
                            return sign * prod(dens[r.reactant_species]) *
                            r.rate_coefficient(temp)
                        end
                    )
                end
            end
        end

        if (s.has_wall_loss)
            push!(sdens_funct_list,
                function (dens::Vector{Float64}, temp::Vector{Float64})
                    return DensWallFluxFunction(dens, temp, s, species_list,
                        reaction_list, system)
                end
            )
        end
    end
    return sdens_funct_list
end


function GenerateTempRateFunctionList(system::System,
    species_list::Vector{Species}, reaction_list::Vector{Reaction})
    # Generate temperature equations
    # Loop over all species
    temp_rate_funct_list = []
    for s in species_list
        stemp_funct_list = GenerateTempRateFunction(s, species_list,
            reaction_list, system)
        if (length(stemp_funct_list)>0)
            push!(temp_rate_funct_list, stemp_funct_list)
        end
    end
    return temp_rate_funct_list
end


function GenerateTempRateFunction(s::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

    stemp_funct_list = Function[]
    if (s.has_temp_eq)
        s_id = s.id
        # "Constants" that are used later
        Q0(dens::Vector{Float64}) = 3.0/2.0 * kb * dens[s_id]
        Q1(temp::Vector{Float64}) = -3.0/2.0 * kb * temp[s_id]

        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s_id
            s_index = findall( x -> x == s_id, r.involved_species )
            if !(length(s_index) == 0)
                # Species s is involved in reaction r
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                if (sign != 0)
                    # Temp gain/loss term is added if there is a non-zero particle balance
                    push!(stemp_funct_list,
                        function (dens::Vector{Float64}, temp::Vector{Float64})
                            return sign * prod(dens[r.reactant_species]) *
                                r.rate_coefficient(temp) *
                                Q1(temp) / Q0(dens)
                        end
                    )
                else
                    # Add energy term due to elastic collisions
                    if (r.id == r_elastic_id)
                        # Find the neutral collision partner
                        neutral_id = r.neutral_species_id 
                        # Set the neutral species mass
                        m_neutral = species_list[neutral_id].mass 
                        m_charged = s.mass
                        Q2 = -3 * kb * m_charged / m_neutral
                        push!(stemp_funct_list,
                            function (dens::Vector{Float64}, temp::Vector{Float64})
                                t_neutral = species_list[r.neutral_species_id].temp0
                                e_rate = Q2 * prod(dens[r.reactant_species]) *
                                    r.rate_coefficient(temp) * (temp[s_id] - t_neutral)
                                return e_rate / Q0(dens)
                            end
                        )
                    end

                end

                # COLLISION INTRINSIC ENERGY GAIN/LOSS
                Er = r.E_threshold
                if (Er != 0) 
                    push!(stemp_funct_list,
                        function (dens::Vector{Float64}, temp::Vector{Float64})
                            e_rate = - Er * prod(dens[r.reactant_species]) *
                                r.rate_coefficient(temp)
                            return e_rate / Q0(dens)
                        end
                    )
                end
            end
        end

        if (s.has_wall_loss)
            push!(stemp_funct_list,
                function (dens::Vector{Float64}, temp::Vector{Float64})
                    flux = TempWallFluxFunction(dens, temp, s, species_list,
                        reaction_list, system)
                    return flux / Q0(dens)
                end
            )
        end

        if (s.has_heating_mechanism)
            push!(stemp_funct_list,
                function (dens::Vector{Float64}, temp::Vector{Float64})
                    ipower = PowerInputFunction(s, system) 
                    return ipower / Q0(dens)
                end
            )
        end
    end
    return stemp_funct_list
end


###############################################################################
# WALL FLUX EXPRESSIONS

function DensWallFluxFunction(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

    if (species.charge == 0)
        dens_wf = 0
    else
        if (system.power_input_method == p_ccp_id)
            wf = WallFluxFunction_CCP(dens, temp, species, species_list,
            reaction_list, system)
            dens_wf = -system.A / system.V * wf
        else
            dens_wf = 0
            print("Flux models only implemented for CCPs\n")
        end
    end
    return dens_wf
end


function TempWallFluxFunction(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction},
    system::System)
    # Input: Species type (must be an ion!)

    if (species.charge == 0)
        temp_wf = 0
    else
        if (system.power_input_method == p_ccp_id)
            if (species.id == s_electron_id)
                # Must account for the electron flux and for the ion flux
                temp_wf = 0
                V_sheath = MeanSheathVoltage(dens, system)
                Te = temp[s_electron_id]
                A = system.A
                V = system.V
                for s in species_list
                    if (s.id == s_electron_id)
                        # Electron flux
                        wf = WallFluxFunction_CCP(dens, temp, s, species_list,
                            reaction_list, system)
                        temp_wf -= wf * A / V * 2.0 * kb * Te
                    elseif (s.charge>0)
                        # Ion flux
                        q = s.charge
                        wf = WallFluxFunction_CCP(dens, temp, s, species_list,
                            reaction_list, system)
                        temp_wf -= wf * A / V * (0.5*kb*Te + q*V_sheath)
                    end
                end
            else
                print("***WARNING*** Temperature flux models for ions are not",
                    " implemented yet\n")
                temp_wf = 0
            end
        else
            print("***WARNING*** Flux models only implemented for CCPs\n")
            temp_wf = 0
        end
    end
    return temp_wf

end


function WallFluxFunction_CCP(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)
    # Input: ions species
    #        if species is electron then summ of all ion species

    if (species.id == s_electron_id)
        # Electron wall flux = sum of all ion wall fluxes
        wallflux = 0
        for s in species_list
            if s.charge>0
                wallflux += s.charge * WallFluxFunction_CCP_ion_species(dens,
                    temp, s, reaction_list, system)
            end
        end
        wallflux = wallflux / abs(species.charge)
    else
        wallflux = WallFluxFunction_CCP_ion_species(dens, temp, species,
            reaction_list, system)
    end
    return wallflux
end


function WallFluxFunction_CCP_ion_species(dens::Vector{Float64},
    temp::Vector{Float64}, s::Species, reaction_list::Vector{Reaction},
    system::System)
    # Input: Species type (must be an ion!)
    
    # Plasma parameters 
    Te = temp[s_electron_id]
    ni = dens[s.id]
    mi = s.mass

    # Ion total collision frequency (elastic collision)
    r_elastic = reaction_list[s.r_elastic_id]
    ng = dens[r_elastic.neutral_species_id]
    K_elastic = r_elastic.rate_coefficient(temp)
    nu = ng * K_elastic

    wallflux = ni * π * kb * Te / system.l / mi / nu
    return wallflux

end


function MeanSheathVoltage(dens::Vector{Float64}, system::System)

    ne = dens[s_electron_id]
    V_sheath = 3/4 * (system.drivI/system.A)^2 / (e * eps0 * ne * system.drivOmega^2)
    return V_sheath
end


end