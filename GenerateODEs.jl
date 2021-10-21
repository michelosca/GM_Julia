module GenerateODEs

using SharedData
using SharedData: A, V, l, drivOmega, drivI
using SharedData: kb, e, eps0
using SharedData: s_electron_id
using SharedData: Species
using SharedData: species_list, reaction_list
using PowerInput: PowerInputFunction

# Symbol flags
global Te_symbol = false
global Te_eV_symbol = false
global m_Ar_symbol = false

# Symbols
global Te_eV = 0.0
global Te = 0.0
global m_Ar = 0.0

function GenerateDensRateFunction(s::Species)

    # This function generates the function terms required in the density ODE function
    # This is for a sigle species s
    # Returns a list of functions. Each of these functions is a term in the ODE
    sdens_funct_list = Function[]

    # Check if species is meant to have a density ODE 
    if (s.has_dens_eq)
        # Loop over the reaction set
        for r in reaction_list
            # Check whether reaction r involves species s 
            s_index = findall( x -> x == s.id, r.involved_species )
            if !(length(s_index) == 0)
                s_index = s_index[1]

                # PARTICLE PRODUCTION rates
                # Terms due to particle gain/loss, e.g. recombination, ionization
                sign = r.species_balance[s_index]
                # Only a density gain/loss term is added if there is a non-zero
                # particle balance
                if !(sign == 0)
                    push!(sdens_funct_list, (dens::Vector{Float64},
                        temp::Vector{Float64}) ->
                        sign * prod(dens[r.reactant_species]) *
                        eval(r.rate_coefficient))
                end
            end
        end

        if (s.has_wall_loss)
            push!(sdens_funct_list, (dens::Vector{Float64},
                temp::Vector{Float64}) -> DensWallFluxFunction(dens, temp, s))
        end
    else
        push!(sdens_funct_list, () -> 0)
    end
    return sdens_funct_list
end


function GenerateDensRateFunctionList()
    # Gathers the Density rate functions from GenerateDensRateFunction for 
    #each species in a list
    
    # Loop over all species
    dens_rate_funct_list = []
    for s in species_list
        sdens_funct_list = GenerateDensRateFunction(s)
        if (length(sdens_funct_list)>0)
            push!(dens_rate_funct_list, sdens_funct_list)
        end
    end
    return dens_rate_funct_list
end


function GenerateTempRateFunction(s::Species)

    stemp_funct_list = Function[]
    if (s.has_temp_eq)
        s_id = s.id
        # "Constants" that are used later
        Q0(dens::Vector{Float64}) = 3.0/2.0 * kb * dens[s_id]
        Q1(temp::Vector{Float64}) = 3.0/2.0 * kb * temp[s_id]

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
                    push!(stemp_funct_list, (dens::Vector{Float64},
                        temp::Vector{Float64}) ->
                        sign * prod(dens[r.reactant_species]) *
                        eval(r.rate_coefficient) * Q1(temp) )
                else
                    # Add energy term due to elastic collisions
                    if (r.id == SharedData.r_elastic_id)
                        # Find the neutral collision partner
                        neutral_id = r.neutral_species_id 
                        # Set the neutral species mass
                        m_neutral = species_list[neutral_id].mass 
                        m_charged = s.mass
                        Q2 = -3 * kb * m_charged / m_neutral
                        push!(stemp_funct_list, (dens::Vector{Float64},
                            temp::Vector{Float64}) -> Q2 *
                            prod(dens[r.reactant_species]) *
                            eval(r.rate_coefficient) *
                            (temp[s_id] - temp[neutral_id]) )
                    end

                end

                # COLLISION INTRINSIC ENERGY GAIN/LOSS
                Er = r.E_threshold
                if (Er != 0) 
                    push!(stemp_funct_list, (dens::Vector{Float64},
                        temp::Vector{Float64}) -> - Er *
                        prod(dens[r.reactant_species]) *
                        eval(r.rate_coefficient) )
                end
            end
        end

        if (s.has_wall_loss)
            push!(stemp_funct_list, (dens::Vector{Float64},
                temp::Vector{Float64}) -> TempWallFluxFunction(dens, temp, s))
        end

        if (s.has_heating_mechanism)
            push!(stemp_funct_list, (dens::Vector{Float64},
                temp::Vector{Float64}) -> PowerInputFunction(dens, temp, s))
        end
    else
        push!(stemp_funct_list, () -> 0)
    end
    return stemp_funct_list
end


function GenerateTempRateFunctionList()
    # Generate temperature equations
    # Loop over all species
    temp_rate_funct_list = []
    for s in species_list
        stemp_funct_list = GenerateTempRateFunction(s)
        if (length(stemp_funct_list)>0)
            push!(temp_rate_funct_list, stemp_funct_list)
        end
    end
    return temp_rate_funct_list
end


function GatherListOfFunctions(dens::Vector{Float64}, temp::Vector{Float64},
    funct_list::Vector{Function})

    f_out = 0
    
    for current_f in funct_list
        f_out += current_f(dens, temp)
    end

    return f_out
end

function UpdateGlobalSymbols(temp::Vector{Float64})
    if Te_eV_symbol
        global Te_eV = temp[SharedData.s_electron_id] * SharedData.K_to_eV
    end

    if Te_symbol
        global Te = temp[SharedData.s_electron_id]
    end
end


function FirstUpdateGlobalSymbols(temp::Vector{Float64})
    if Te_eV_symbol
        global Te_eV = temp[SharedData.s_electron_id] * SharedData.K_to_eV
        print("Te_eV ", Te_eV, "\n")
    end

    if Te_symbol
        global Te = temp[SharedData.s_electron_id]
    end

    if m_Ar_symbol
        global m_Ar = 4 * SharedData.mp
    end
end

function SetSymbolFlag(num::Int64)

    errcode = 0

    if (num == 1)
        global Te_symbol = true
    elseif (num == 2)
        global Te_eV_symbol = true
    elseif (num == 3)
        global m_Ar_symbol = true
    else
        print("***ERROR*** The parameter in rate coefficient is not recognized\n")
        errcode = 1
    end

    return errcode
end


###############################################################################
# WALL FLUX EXPRESSIONS

function WallFluxFunction_CCP(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species)
    # Input: ions species
    #        if species is electron then summ of all ion species

    if (species.id == s_electron_id)
        # Electron wall flux = sum of all ion wall fluxes
        wallflux = 0
        for s in SharedData.species_list
            if s.charge>0
                wallflux += s.charge * WallFluxFunction_CCP_ion_species(dens, temp, s)
            end
        end
        wallflux = wallflux / abs(species.charge)
    else
        wallflux = WallFluxFunction_CCP_ion_species(dens, temp, species)
    end
    return wallflux
end


function WallFluxFunction_CCP_ion_species(dens::Vector{Float64},
    temp::Vector{Float64}, s::Species)
    # Input: Species type (must be an ion!)
    
    # Plasma parameters 
    Te = temp[s_electron_id]
    ni = dens[s.id]
    mi = s.mass

    # Ion total collision frequency (elastic collision)
    r_elastic = SharedData.reaction_list[s.r_elastic_id]
    ng = dens[r_elastic.neutral_species_id]
    K_elastic = eval(r_elastic.rate_coefficient)
    nu = ng * K_elastic

    wallflux = ni * π * kb * Te / l / mi / nu
    return wallflux

end


function MeanSheathVoltage(dens::Vector{Float64})

    ne = dens[s_electron_id]
    V_sheath = 3/4 * (drivI/A)^2 / (e * eps0 * ne * drivOmega^2)
    return V_sheath
end


function TempWallFluxFunction(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species)
    # Input: Species type (must be an ion!)

    if (species.charge == 0)
        temp_wf = 0
    else
        if (SharedData.power_input_method == SharedData.p_ccp_id)
            if (species.id == s_electron_id)
                # Must account for the electron flux and for the ion flux
                temp_wf = 0
                V_sheath = MeanSheathVoltage(dens)
                Te = temp[s_electron_id]
                for s in SharedData.species_list
                    if (s.id == s_electron_id)
                        # Electron flux
                        wf = WallFluxFunction_CCP(dens, temp, s)
                        Te = temp[s_electron_id]
                        temp_wf -= wf * A / V * 2 * kb * Te
                    elseif (s.charge>0)
                        # Ion flux
                        q = s.charge
                        wf = WallFluxFunction_CCP(dens, temp, s)
                        temp_wf -= wf * A / V * (0.5*kb*Te + q*V_sheath)
                    end
                end
            else
                print("Temperature flux models for ions have not been implemented\n")
                temp_wf = 0
            end
        else
            print("Flux models only implemented for CCPs\n")
            temp_wf = 0
        end
    end
    return temp_wf

end

function DensWallFluxFunction(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species)

    if (species.charge == 0)
        dens_wf = 0
    else
        if (SharedData.power_input_method == SharedData.p_ccp_id)
            wf = WallFluxFunction_CCP(dens, temp, species)
            dens_wf = -A / V * wf
        else
            dens_wf = 0
            print("Flux models only implemented for CCPs\n")
        end
    end
    return dens_wf
end

end