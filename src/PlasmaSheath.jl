# Copyright (C) 2021 Michel Osca Engelbrecht
#
# This file is part of GM Julia.
#
# GM Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GM Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module PlasmaSheath

using PowerInput: PowerInputFunction
using SharedData: Species, Reaction, System, SpeciesID
using SharedData: me, K_to_eV, amu, e, eps0
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation
using Roots: find_zeros

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GetSheathVoltage
#   - SheathVoltage_FluxBalanceEquation
#   - SheathVoltage_OhmicPowerSource
#   - SheathVoltage_InterpolateFluxEquation 

function GetSheathVoltage!(system::System, species_list::Vector{Species},
    sID::SpeciesID, time::Float64)

    errcode = 0

    # Obtain the positive charged particles flux and solve for the potential
    solve_method = system.Vsheath_solving_method

    V_sheath = 0.0
    if solve_method == s_ohmic_power 
        electron_species = species_list[sID.electron]
        V_sheath = SheathVoltage_OhmicPowerSource(electron_species, system,
            sID, time)
    elseif solve_method == s_flux_balance 
        V_sheath = SheathVoltage_FluxBalanceEquation(species_list,
            sID.electron)
    elseif solve_method == s_flux_interpolation
        V_sheath = SheathVoltage_InterpolateFluxEquation(species_list, system)
    else
        print("***WARNING*** No potential sheath calculation done\n")
    end
    system.plasma_potential = V_sheath
    return errcode
end

function SheathVoltage_FluxBalanceEquation(species_list::Vector{Species},
    electron_id::Int64)

    positive_flux = 0.0 # Fluxes of positive charged particles
    electrons = species_list[electron_id]
    n_e = electrons.dens
    Te_eV = electrons.temp * K_to_eV
    v_th = electrons.v_thermal
    for s in species_list
        if s.charge > 0.0
            positive_flux += s.flux 
        end
    end
    #print("Positive flux ", positive_flux,"\n")
    V_sheath = -log(positive_flux / v_th / n_e * 4.0) * Te_eV

    return V_sheath
end

function SheathVoltage_OhmicPowerSource(electrons::Species, system::System, sID::SpeciesID, time::Float64)
    # From: A. Hurlbatt et al. (2017)

    S_ohm = PowerInputFunction(electrons, system, sID, time)
    mfp = electrons.mfp 
    v_th = electrons.v_thermal
    V_sheath = 3.0/2.0 * S_ohm*e*mfp / (me*eps0*system.drivOmega^2*v_th)

    return V_sheath
end


function SheathVoltage_InterpolateFluxEquation(species_list::Vector{Species},
    system::System)

    positive_flux = 0.0        # Fluxes of positive charged particles
    for s in species_list
        if s.charge > 0.0
            positive_flux += s.flux
        end
    end

    function negative_flux_funct(pot::Float64)
        negative_flux = 0.0        # Fluxes of negative charged particles
        for s in species_list
            if s.charge < 0.0
                v_th = s.v_thermal
                T_eV = s.temp * K_to_eV
                negative_flux += 0.25 * s.dens * v_th * exp(-pot/T_eV) 
            end
        end
        return negative_flux - positive_flux
    end

    # Solve for plasma potential 
    interp_flag = true
    fact = 1.5
    V_guess = copy(system.plasma_potential)
    while interp_flag
        V_guess *= fact 
        V_sheath = find_zeros(negative_flux_funct, (0,V_guess))
        if length(V_sheath) == 1
            interp_flag = false
            V_guess = V_sheath[1]
        elseif length(V_sheath) < 1 
            fact = 1.5
        elseif length(V_sheath) > 1 
            fact = 0.5
        end
    end
    return V_guess
end

end