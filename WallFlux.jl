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

module WallFlux

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: kb, e, K_to_eV 
using SharedData: p_ccp_id, p_icp_id 
using SharedData: s_ohmic_power

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - DensWallFluxFunction
# - TempWallFluxFunction
# - GetIonFlux
#   - IonParticleFlux_ICP
#   - IonParticleFlux_CCP
# - GetElectronFlux!

function DensWallFluxFunction(species::Species, system::System)
    A = system.A
    V = system.V
    particle_flux = -A /V * species.flux
    return particle_flux 
end

function TempWallFluxFunction(temp::Vector{Float64}, species::Species,
    species_list::Vector{Species}, system::System, V_sheath::Float64,
    sID::SpeciesID)

    A = system.A
    V = system.V
    s_id = species.id
    temp_flux = 0.0

    if s_id == sID.electron
        Te = temp[s_id]
        temp_flux -= 2.0 * kb * Te * species.flux
        # Add ion flux terms
        for s in species_list
            charge = s.charge
            if (charge > 0)
                temp_flux -= (0.5*kb*Te + charge*V_sheath) * s.flux
            end
        end
    end

    temp_flux *= A/V*2
    return temp_flux
end


function UpdatePositiveFlux!(species_list::Vector{Species}, system::System)

    # Geometrical factors
    if system.power_input_method == p_icp_id
        L = system.l
        R = system.radius
    end

    # Particle flux of positive ions
    for s in species_list
        if s.charge > 0
            if system.power_input_method == p_ccp_id
                s.flux = s.v_Bohm * s.n_sheath 
                #print("Spcies ", s.id, " flux ", s.flux,"\n")
            elseif system.power_input_method == p_icp_id
                h_L = s.h_L
                h_R = s.h_R
                fac = (R^2 * h_L + R*L*h_R) / (R^2 + R*L)
                s.flux = fac * s.v_Bohm * s.n_sheath
            end
        end
    end
end


function UpdateNegativeFlux!(species_list::Vector{Species}, system::System,
    sID::SpeciesID, V_sheath::Float64)

    # Electron flux solved from the flux balance equation, i.e.
    # flux_electrons = SUM(flux_ions)
    # Negative ion species are assumed to have zero net flow through the walls

    solve_method = system.Vsheath_solving_method
    if solve_method == s_ohmic_power 
        negative_flux = 0.0
    end

    for s in species_list
        if solve_method == s_ohmic_power 
            q = s.charge
            if q > 0
                negative_flux += s.flux * q 
            end
        else
            if s.charge < 0
                T_eV = s.temp * K_to_eV
                s.flux = 0.25 * s.n_sheath * s.v_thermal * exp(-V_sheath / T_eV)
            end
        end
    end
    if solve_method == s_ohmic_power 
        species_list[sID.electron].flux = negative_flux / e 
    end

end

end