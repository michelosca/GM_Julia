module PlasmaSheath

using PowerInput: PowerInputFunction
using SharedData: Species, Reaction, System, SpeciesID
using SharedData: me, K_to_eV, amu, e, eps0
using SharedData: p_icp_id, p_ccp_id
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation
using Roots: find_zeros

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GetSheathVoltage
#   - SheathVoltage_FluxBalanceEquation
#   - SheathVoltage_OhmicPowerSource
#   - SheathVoltage_InterpolateFluxEquation 

function GetSheathVoltage(species_list::Vector{Species}, system::System, sID::SpeciesID)

    # Obtain the positive and negatice charged particles flux and solve for the potential
    solve_method = system.Vsheath_solving_method

    if solve_method == s_ohmic_power 
        electron_species = species_list[sID.electron]
        V_sheath = SheathVoltage_OhmicPowerSource(electron_species, system, sID)
    elseif solve_method == s_flux_balance 
        V_sheath = SheathVoltage_FluxBalanceEquation(species_list, sID.electron)
    elseif solve_method == s_flux_interpolation
        V_sheath = SheathVoltage_InterpolateFluxEquation(species_list)
    else
        print("***WARNING*** No potential sheath calculation done\n")
        V_sheath = 0.0 
    end
    return V_sheath
end

function SheathVoltage_FluxBalanceEquation(species_list::Vector{Species},
    electron_id::Int64)

    positive_flux = 0.0 # Fluxes of positive charged particles
    electrons = species_list[electron_id]
    n_e = electrons.n_sheath 
    Te_eV = electrons.temp * K_to_eV
    v_th = electrons.v_thermal
    for s in species_list
        if s.charge > 0
            positive_flux += s.flux 
        end
    end
    #print("Positive flux ", positive_flux,"\n")
    V_sheath = -log(positive_flux / v_th / n_e * 4.0) * Te_eV

    return V_sheath
end

function SheathVoltage_OhmicPowerSource(electrons::Species, system::System, sID::SpeciesID)
    # From: A. Hurlbatt et al. (2017)

    S_ohm = PowerInputFunction(electrons, system, sID)
    mfp = electrons.mfp 
    v_th = electrons.v_thermal
    V_sheath = 3.0/2.0 * S_ohm*e*mfp / (me*eps0*system.drivOmega^2*v_th)

    return V_sheath
end


function SheathVoltage_InterpolateFluxEquation(species_list::Vector{Species})

    positive_flux = 0.0        # Fluxes of positive charged particles
    for s in species_list
        if s.charge > 0
            positive_flux += s.flux
        end
    end

    function negative_flux_funct(pot::Float64)
        negative_flux = 0.0        # Fluxes of negative charged particles
        for s in species_list
            if s.charge < 0
                n_sheath = s.n_sheath 
                v_th = s.v_thermal
                T_eV = s.temp * K_to_eV
                negative_flux += 0.25 * n_sheath * v_th * exp(-pot/T_eV) 
            end
        end
        return negative_flux - positive_flux
    end

    # Solve now plasma potential 
    V_sheath = find_zeros(negative_flux_funct, (-1000,1000))
    if (length(V_sheath) > 1 || V_sheath == Float64[])
        print("***ERROR*** V_sheath interpolation problem. Length ", length(V_sheath),"\n")
    end
    return V_sheath[1]
end

end