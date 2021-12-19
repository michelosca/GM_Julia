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

module PlasmaParameters

using SharedData: Species, Reaction, System, SpeciesID
using SharedData: kb, K_to_eV, e
using SharedData: p_icp_id, p_ccp_id
using EvaluateExpressions: ReplaceExpressionValues

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - UpdateSpeciesParameters
# - GetBohmSpeed
# - GetThermalSpeed
# - GetMFP

function UpdateSpeciesParameters!(temp::Vector{Float64}, dens::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)

    # Update first dens and temperature
    for s in species_list
        s.temp = temp[s.id]
        s.dens = dens[s.id]
    end

    for s in species_list
        s.v_thermal = GetThermalSpeed(s)
        s.v_Bohm = GetBohmSpeed(temp[sID.electron], s.mass)
        s.mfp = GetMFP(temp, dens, s, species_list, system, sID)

        if system.power_input_method == p_ccp_id
            if s.charge > 0
                s.n_sheath = GetSheathDensity(s, system)
            else
                s.n_sheath = s.dens
            end
        elseif system.power_input_method == p_icp_id
            UpdateAlpha!(dens, species_list, system, sID.electron)

            s.n_sheath = s.dens
            if (s.id == sID.electron)
                continue
            end

            if s.charge > 0
                s.h_L, s.h_R = Get_h_Parameters(temp, dens, s, species_list,
                    system, sID.electron)
            end

            s.D = GetNeutralDiffusionCoeff(s, species_list, sID)
        end
    end
end


function GetMFP(temp::Vector{Float64}, dens::Vector{Float64}, species::Species,
    species_list::Vector{Species}, system::System, sID::SpeciesID)

    ilambda = 0.0         # inverse mean-free-path

    v_th_s = species.v_thermal
    id = species.id

    for r in species.reaction_list

        # Get max. thermal speed of reacting species
        v_th = v_th_s
        for i in r.reactant_species
            v_th_n = species_list[i].v_thermal
            v_th = max(v_th_n, v_th)
        end
        
        # Exclude the species for which the mfp is calculated
        r_species = copy(r.reactant_species)
        index = findall( x -> x == id, r_species )
        deleteat!(r_species, index)
        
        # Set collision cross section
        if system.prerun
            K = r.rate_coefficient(temp, sID)
        else
            K = ReplaceExpressionValues(r.rate_coefficient, temp,
                species_list, system, sID)
        end
        cross_section = K / v_th

        # Density of colliding partners
        n = prod(dens[r_species])

        # Add to mfp-buffer
        ilambda += n * cross_section 
    end

    if ilambda == 0
        ilambda = 1.e-100
    end

    lambda = 1.0/ilambda
    return lambda
end


function UpdateAlpha!(dens::Vector{Float64}, species_list::Vector{Species},
    system::System, electron_id::Int64)

    # Alpha parameter
    alpha = 0.0
    for s in species_list
        q = s.charge
        if s.id == electron_id
            continue
        elseif q < 0
            alpha += dens[s.id]
        end
    end
    alpha /= dens[electron_id]
    system.alpha = alpha
end


function GetBohmSpeed(Te::Float64, mass::Float64)

    uB = sqrt(kb * Te / mass)
    return uB
end


function GetThermalSpeed(species::Species)

    T = species.temp
    m = species.mass
    v_th = sqrt(8.0* kb * T /  m / pi)
    return v_th
end


function Get_h_Parameters(temp::Vector{Float64}, dens::Vector{Float64}, species::Species,
    species_list::Vector{Species}, system::System, electron_id::Int64)
    # See Gudmundsson (2000) On the plasma parameters of a planar inductive oxygen discharge
    # Only valid for positive charged species

    L = system.l
    R = system.radius

    # Gamma parameter
    Te = temp[electron_id]
    Ti = temp[species.id]
    gamma = Te / Ti
    alpha = system.alpha

    # Get species mean free path
    lambda = species.mfp
    h_L = 0.86 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(3.0+0.5*L/lambda)
    h_R = 0.8 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(4.0+R/lambda)
    return h_L, h_R
end


function GetLambda(system::System)

    L = system.l
    R = system.radius
    Lambda = 1.0/sqrt((pi/L)^2 + (2.405/R)^2)

    return Lambda
end


function GetNeutralDiffusionCoeff(species::Species,
    species_list::Vector{Species}, sID::SpeciesID)
    # Neutral Diffusion Coefficient

    T_eV = species.temp * K_to_eV
    mfp = species.mfp
    D = e * T_eV * mfp / (species.v_thermal * species.mass)
    return D
end

function GetSheathDensity(species::Species, system::System)

    ni = species.dens
    uB = species.v_Bohm 
    uTh = species.v_thermal 
    mfp = species.mfp
    l = system.l # System length

    n_sheath = pi * ni * (uB / uTh) * (mfp / l) # Density at sheath edge
    return n_sheath
end

end