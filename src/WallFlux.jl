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
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation

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
    particle_flux = -A/V * species.flux
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
        temp_flux -= 2.0*kb*Te * species.flux
        # Add ion flux terms
        for s in species_list
            charge = s.charge
            if charge > 0.0
                temp_flux -= (0.5*kb*Te + charge*V_sheath) * s.flux
            end
        end
    end

    temp_flux *= A/V
    return temp_flux
end


function UpdatePositiveFlux!(species_list::Vector{Species})

    errcode = 0

    # Particle flux of positive ions
    for s in species_list
        if s.charge > 0.0
            flux = s.v_Bohm * s.n_sheath
            s.flux = flux
            # Charged species flux recombines, and it's corresponding
            # neutral species gains this mass
            species_list[s.species_id].flux += -flux
        end
    end
    return errcode
end


function UpdateNegativeFlux!(species_list::Vector{Species}, system::System,
    sID::SpeciesID)

    errcode = 0

    solve_method = system.Vsheath_solving_method
    if solve_method == s_ohmic_power 
        # This model assumes there is no mass flow of neg-Ions through the
        #sheath to the wall, i.e. flux_electrons = SUM(flux_ions)
        negative_flux = 0.0
    end

    if solve_method == s_flux_balance
        # This methods assumes there is no flux from negative ions
        electrons = species_list[sID.electron]
        T_eV = electrons.temp * K_to_eV
        n_e = electrons.dens
        electrons.flux = 0.25 * n_e * electrons.v_thermal *
            exp(-system.plasma_potential/ T_eV)
    else
        for s in species_list
            q = s.charge
            if solve_method == s_ohmic_power 
                if q > 0.0
                    negative_flux += s.flux * q 
                end
            elseif solve_method == s_flux_interpolation
                if q < 0.0
                    T_eV = s.temp * K_to_eV
                    s.flux = 0.25 * s.dens * s.v_thermal *
                        exp(-system.plasma_potential/ T_eV)
                end
            end
        end

        if solve_method == s_ohmic_power 
            species_list[sID.electron].flux = negative_flux / e 
        end
    end

    return errcode
end

end