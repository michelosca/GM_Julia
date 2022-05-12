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
using SharedData: h_classical, h_Gudmundsson, h_Monahan 
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

    # FIRST: Update dens and temperature
    for s in species_list
        s.temp = temp[s.id]
        s.dens = dens[s.id]
        s.pressure = s.dens * kb * s.temp
    end

    # SECOND: Update parameters that depend on dens,temp and other species parameters
    system.total_pressure = UpdateTotalPressure(species_list, sID)
    if system.h_id == h_Gudmundsson || system.h_id == h_Monahan 
        system.alpha = UpdateAlpha(dens, species_list, sID.electron)
    end

    for s in species_list
        s.v_thermal = GetThermalSpeed(s)
        s.v_Bohm = GetBohmSpeed(temp[sID.electron], s.mass)
        s.mfp = GetMFP(temp, dens, s, species_list, system, sID)

        s.gamma = GetStickingCoefficient(s, species_list, sID)
        s.D = GetNeutralDiffusionCoeff(s)

        if s.charge > 0.0
            s.n_sheath = GetSheathDensity(s, temp, species_list, system, sID)
        end
    end
end

function UpdateTotalPressure(species_list::Vector{Species}, sID::SpeciesID)
    p_total = 0.0
    for s in species_list
        if s.id == sID.electron
            continue
        end
        p_total += s.pressure 
    end
    return p_total
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
        if n < 0.0
            n = 0.0
        end

        # Add to mfp-buffer
        ilambda += n * cross_section 
    end

    if ilambda == 0
        ilambda = 1.e-100
    end

    lambda = 1.0/ilambda
    return lambda
end


function UpdateAlpha(dens::Vector{Float64}, species_list::Vector{Species},
    electron_id::Int64)
    # Alpha parameter: is the ratio of negative ion species to electron density

    alpha = 0.0
    for s in species_list
        if s.id == electron_id
            continue
        elseif s.charge < 0.0
            n = dens[s.id]
            if n > 0.0
                alpha += n 
            end
        end
    end
    alpha /= dens[electron_id]
    return alpha
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


function GetLambda(system::System)

    L = system.l
    R = system.radius
    Lambda = 1.0/sqrt((pi/L)^2 + (2.405/R)^2)

    return Lambda
end


function GetNeutralDiffusionCoeff(species::Species)
    # Neutral Diffusion Coefficient

    mfp = species.mfp
    D = e * species.temp * mfp / (species.v_thermal * species.mass)
    return D
end

function GetSheathDensity(species::Species, temp::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)
    # Calculates the density at sheath edge
    # This is only interesting for positive ions. Used later on to
    # calculate the wall-flux = u_B * n_Sheath

    n_0 = species.dens

    if system.h_id == h_classical
        # Electropositive plasmas
        L = system.l
        mfp = species.mfp

        #if system.total_pressure < 0.26664474 # ( = 2 mTorr) 
        #    # Low pressure range
        #    h = 0.425 
        #elseif system.total_pressure < 20 # ( = 150 mTorr)
        #    # Intermediate pressure range
        #    h = 0.86 / sqrt(3.0 + L/mfp)
        #else
           # High pressure
           uB = species.v_Bohm 
           uTh = species.v_thermal 
           h = pi * (uB / uTh) * (mfp / L)
        #end
    elseif system.h_id == h_Gudmundsson
        # ICP Ar/O2
        R = system.radius
        L = system.l
        Te = temp[sID.electron]
        Ti = temp[species.id]
        gamma = Te / Ti
        alpha = system.alpha
        lambda = species.mfp

        h_L = 0.86 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(3.0+0.5*L/lambda)
        h_R = 0.8 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(4.0+0.5*R/lambda)

        h = (R^2 * h_L + R*L*h_R) / (R^2 + R*L)

    elseif system.h_id == h_Monahan
        #print("NAME: ", species.name,"\n")
        #print("n_0 " ,n_0,"\n")
        if n_0 > 0.0
            # Electronegative plasmas
            L = system.l
            alpha = system.alpha
            Te = species_list[sID.electron].temp
            Ti = species.temp
            sqrt_Te_Ti = sqrt(Te/Ti)
            mfp = species.mfp
            uTh = species.v_thermal 
            K_recombination = GetRecombinationRate(species, temp, species_list, system, sID)
            ni_star = 15.0/56.0 * uTh / K_recombination / mfp
            n_n0_p3_2 = (alpha * species_list[sID.electron].dens)^1.5

            h_a = 0.86 / sqrt(3.0 + L/mfp) / (1.0 + alpha) 
            h_b = alpha / (1.0 + alpha) / ( sqrt_Te_Ti * (1.0+1.0/sqrt(2.0*pi)/mfp) )
            h_c = 1.0/( sqrt_Te_Ti * (1.0+sqrt(ni_star)*n_0/n_n0_p3_2) ) 
            h = sqrt(h_a^2 + h_b^2 + h_c^2)

        #    print("sqrt_Te_Ti " ,sqrt_Te_Ti,"\n")
        #    print("mfp " ,mfp ,"\n")
        #    print("K_rec " ,K_recombination ,"\n")
        #    print("ni_star " ,ni_star ,"\n")
        #    print("n_n0_p3_2 " ,n_n0_p3_2 ,"\n")
        #    print("h_a " , h_a,"\n")
        #    print("h_b " , h_b,"\n")
        #    print("h_c " , h_c,"\n")
        #    print("h " , h,"\n\n")
        else
            h = 0.0
        end
    end
    n_sheath = n_0 * h

    return n_sheath
end


function GetRecombinationRate(species::Species, temp::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)

    s_id = species.id

    K_recombination = 0.0
    for r in species.reaction_list
        s_index = findall( x -> x == s_id, r.involved_species )[1]
        sign = r.species_balance[s_index]
        if sign < 0
            if system.prerun
                K = r.rate_coefficient(temp, sID) 
            else
                K = ReplaceExpressionValues(r.rate_coefficient, temp,
                    species_list, system, sID)
            end
            K_recombination += K
        end
    end
    return K_recombination
end


function GetStickingCoefficient(species::Species,
    species_list::Vector{Species}, sID::SpeciesID)

    gamma = species.gamma

    if species.id == sID.O
        # This is for stainless steel walls ( Gudmundsson 2007)
        s = species_list[sID.O2]
        pO2_mTorr = s.dens * kb * s.temp * 0.13332237

        if pO2_mTorr < 2
            gamma = 1.0 - pO2_mTorr * 0.25
        else
            gamma = 0.1438 * exp(2.5069/pO2_mTorr)
        end
        #print("Sticking coeff. ", gamma,"\n")
    end

    return gamma
end


end