module PlasmaSheath

using PlasmaParameters: GetMFP, GetThermalSpeed, GetBohmSpeed
using PowerInput: PowerInputFunction
using SharedData: Species, Reaction, System
using SharedData: me, K_to_eV, amu, e, eps0
using SharedData: p_icp_id, p_ccp_id
using InputBlock_Species: s_electron_id
using Roots: find_zeros

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GetSheathVoltage
#   - SheathVoltage_FluxBalanceEquation
#   - SheathVoltage_OhmicPowerSource
#   - SheathVoltage_InterpolateFluxEquation 

function GetSheathVoltage(dens::Vector{Float64}, temp::Vector{Float64},
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, flux_list::Vector{Tuple})

    # Obtain the positive and negatice charged particles flux and solve for the potential
    ipower_method = system.power_input_method

    if ipower_method == p_ccp_id
        V_sheath = SheathVoltage_OhmicPowerSource(dens, temp, species_list,
            reaction_list, system)
        #V_sheath = SheathVoltage_FluxBalanceEquation(temp, species_list,
        #    flux_list)
    elseif ipower_method == p_icp_id
        V_sheath = SheathVoltage_FluxBalanceEquation(temp, species_list,
            flux_list)
    else
        print("***WARNING*** No potential sheath calculation done\n")
        V_sheath = 0.0 
    end
    return V_sheath
end

function SheathVoltage_FluxBalanceEquation(temp::Vector{Float64},
    species_list::Vector{Species}, flux_list::Vector{Tuple})

    flux_plus = 0.0 # Fluxes of positive charged particles
    n_sheath_e = 0.0
    for s in species_list
        if s.charge > 0
            flux_plus += flux_list[s.id][1] 
            n_sheath_e += flux_list[s.id][2]
        end
    end
    v_th = GetThermalSpeed(temp, species_list[s_electron_id])
    Te_eV = temp[s_electron_id] * K_to_eV
    V_sheath = -log(flux_plus / v_th / n_sheath_e * 4) * Te_eV

    return V_sheath
end

function SheathVoltage_OhmicPowerSource(dens::Vector{Float64},
    temp::Vector{Float64}, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)
    # From: A. Hurlbatt et al. (2017)

    electrons = species_list[s_electron_id]
    S_ohm = PowerInputFunction(electrons, system)
    e_mfp = GetMFP(dens, temp, electrons, species_list, reaction_list)
    v_th = GetThermalSpeed(temp, electrons)
    V_sheath = 3/2*S_ohm*e*e_mfp/me/eps0/system.drivOmega^2/v_th

    return V_sheath
end


function SheathVoltage_InterpolateFluxEquation( 
    temp::Vector{Float64}, species_list::Vector{Species})

    flux_plus = 0.0        # Fluxes of positive charged particles
    flux_min = Function[]  # Fluxes of netive charged particles
    for s in species_list
        if s.charge > 0
            flux_plus += flux_list[s.id][1]
        elseif s.charge < 0
            push!(flux_min,
                function (pot::Float64)
                    n_sheath = flux_list[s.id][2]
                    v_th = GetThermalSpeed(temp, s) 
                    T_eV = temp[s.id] * K_to_eV
                    V_sheath = 0.25 * n_sheath * v_th * exp(-pot/T_eV) 
                end
            )
        end
    end

    # Solve now "flux_min - flux_plus = 0" and get the plasma potential
    function flux_min_funct(pot::Float64)
        output = 0.0
        for f in flux_min
            output += f(pot)
        end
        return output - flux_plus
    end
    V_sheath = find_zeros(flux_min_funct, (-1000,1000))
    if (length(V_sheath) > 1)
        print("***WARNING*** More than one root found while V_sheath interpolation!\n")
    end
    return V_sheath[1]
end

end