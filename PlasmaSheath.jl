module PlasmaSheath

using PlasmaParameters: GetMFP, GetThermalSpeed, GetBohmSpeed
using PowerInput: PowerInputFunction
using SharedData: Species, Reaction, System
using SharedData: me, K_to_eV, amu, e, eps0
using SharedData: p_icp_id, p_ccp_id
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation
using InputBlock_Species: s_electron_id
using Roots: find_zeros

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GetSheathVoltage
#   - SheathVoltage_FluxBalanceEquation
#   - SheathVoltage_OhmicPowerSource
#   - SheathVoltage_InterpolateFluxEquation 

function GetSheathVoltage(species_list::Vector{Species}, system::System)

    # Obtain the positive and negatice charged particles flux and solve for the potential
    solve_method = system.Vsheath_solving_method

    if solve_method == s_ohmic_power 
        electron_species = species_list[s_electron_id]
        V_sheath = SheathVoltage_OhmicPowerSource(electron_species, system)
    elseif solve_method == s_flux_balance 
        V_sheath = SheathVoltage_FluxBalanceEquation(species_list)
    elseif solve_method == s_flux_interpolation
        V_sheath = SheathVoltage_InterpolateFluxEquation(species_list)
    else
        print("***WARNING*** No potential sheath calculation done\n")
        V_sheath = 0.0 
    end
    return V_sheath
end

function SheathVoltage_FluxBalanceEquation(species_list::Vector{Species})

    flux_plus = 0.0 # Fluxes of positive charged particles
    n_e = 0.0
    Te_eV = 0.0
    v_th = 0.0
    for s in species_list
        if s.charge > 0
            flux_plus += s.flux 
            n_e += s.n_sheath 
        elseif s.id == s_electron_id
            Te_eV += s.temp * K_to_eV
            v_th += s.v_thermal
        end
    end
    V_sheath = -log(flux_plus / v_th / n_e * 4) * Te_eV

    return V_sheath
end

function SheathVoltage_OhmicPowerSource(electrons::Species, system::System)
    # From: A. Hurlbatt et al. (2017)

    S_ohm = PowerInputFunction(electrons, system)
    mfp = electrons.mfp 
    v_th = electrons.v_thermal
    V_sheath = 3.0/2.0 * S_ohm*e*mfp / (me*eps0*system.drivOmega^2*v_th)

    return V_sheath
end


function SheathVoltage_InterpolateFluxEquation(species_list::Vector{Species})

    flux_plus = 0.0        # Fluxes of positive charged particles
    flux_min = Function[]  # Fluxes of netive charged particles
    for s in species_list
        if s.charge > 0
            flux_plus += s.flux
        elseif s.charge < 0
            n_sheath = s.n_sheath 
            v_th = s.v_thermal
            T_eV = s.temp * K_to_eV
            push!(flux_min,
                function (pot::Float64)
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