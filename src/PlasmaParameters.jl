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
using SharedData: kb, K_to_eV, e, r_extended
using SharedData: c_io_error, r_diffusion, r_lower_threshold, r_emission_rate 
using SharedData: h_classical, h_Gudmundsson, h_Monahan 
using EvaluateExpressions: ReplaceExpressionValues
using Printf: @sprintf
using PrintModule: PrintErrorMessage, PrintWarningMessage 

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################

function UpdateParameters!(temp::Vector{Float64}, dens::Vector{Float64},
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, sID::SpeciesID)

    errcode = 0

    # FIRST: Update dens, temperature and pressure. Set flux values to zero.
    for s in species_list

        #print("ne: ", species_list[sID.electron].dens,"; ")
        #print("Te: ", species_list[sID.electron].temp*K_to_eV,"\n")
        # Temperature
        s.temp = temp[s.id]
        if s.temp < 0.0 || isnan(s.temp)
            err_message = @sprintf("%s temperature is negative: %15g eV",
                s.name, s.temp*K_to_eV)
            PrintErrorMessage(system, err_message) 
            return c_io_error
        end

        # Density
        s.dens = dens[s.id]
        if s.id == sID.electron
            if s.dens < 0.0 || isnan(s.dens)
                err_message = @sprintf("%s density is negative: %15g m^-3",
                    s.name, s.dens)
                PrintWarningMessage(system, err_message) 
            end
        end

        # Pressure
        s.pressure = s.dens * kb * s.temp

        # Reset species wall flux
        s.flux = 0.0
    end

    # SECOND: Update rate coefficient values that only depend on temperature
    errcode = UpdateRateCoefficientValues!(reaction_list, dens, temp, species_list,
        system, sID, false)
    if errcode == c_io_error
        PrintErrorMessage(system, "UpdateRateCoefficientValues (regular collisions) failed")
        return c_io_error
    end

    # THIRD: Update parameters that depend on dens, temp and other species parameters
    errcode = UpdateTotalPressure!(system, species_list, sID)
    if errcode == c_io_error
        PrintErrorMessage(system, "UpdateTotalPressure failed")
        return c_io_error
    end

    # FOURTH: Update species parameters
    errcode = UpdateSpeciesParameters!(temp, dens, species_list, system, sID)
    if errcode == c_io_error
        PrintErrorMessage(system, "UpdateSpeciesParameters! failed")
        return c_io_error
    end

    # FIFTH: Update special rate coefficients
    errcode = UpdateRateCoefficientValues!(reaction_list, dens, temp, species_list,
        system, sID, true)
    if errcode == c_io_error
        PrintErrorMessage(system, "UpdateRateCoefficientValues (species collisions) failed")
        return c_io_error
    end
    return errcode 
end


function UpdateSpeciesParameters!(temp::Vector{Float64}, dens::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)
    
    errcode = 0

    if system.h_id == h_Gudmundsson || system.h_id == h_Monahan 
        errcode = UpdateElectronegativity!(system, dens, species_list,
            sID.electron)
        if errcode == c_io_error
            PrintErrorMessage(system, "UpdateElectronegativity failed")
            return c_io_error
        end
    end


    for s in species_list
        errcode = GetThermalSpeed!(s)
        if errcode == c_io_error
            err_message = @sprintf("GetThermalSpeed for %s failed", s.name)
            PrintErrorMessage(system, err_message)
            return c_io_error
        end

        errcode = GetBohmSpeed!(s, temp[sID.electron])
        if errcode == c_io_error
            err_message = @sprintf("GetBohmSpeed for %s failed", s.name)
            PrintErrorMessage(system, err_message)
            return c_io_error
        end
        
        errcode = GetMFP!(s, dens, species_list)
        if errcode == c_io_error
            err_message = @sprintf("GetMFP for %s failed", s.name)
            PrintErrorMessage(system, err_message)
            return c_io_error
        end

        errcode = GetStickingCoefficient!(s, species_list, sID)
        errcode = GetNeutralDiffusionCoeff!(s)
        if errcode == c_io_error
            err_message = @sprintf("GetNeutralDiffusionCoeff for %s failed",
                s.name)
            PrintErrorMessage(system, err_message)
            return c_io_error
        end

        if s.charge > 0.0
            errcode = GetSheathDensity!(s, species_list, system, sID)
            if errcode == c_io_error
                err_message = @sprintf("GetSheathDensity for %s failed",s.name)
                PrintErrorMessage(system, err_message)
                return c_io_error
            end
        end
    end

    return errcode
end


function UpdateRateCoefficientValues!(reaction_list::Vector{Reaction},
    dens::Vector{Float64}, temp::Vector{Float64}, species_list::Vector{Species},
    system::System, sID::SpeciesID, special_coll::Bool)

    errcode = 0

    for r in reaction_list
        # Updates wall-loss reactions
        if r.case == r_diffusion || r.case == r_emission_rate
            if special_coll 
                r.K_value = r.rate_coefficient(dens, temp, species_list,
                    system, sID) 
                
                # Self-absorption correction in emission reactions
                if r.case == r_emission_rate && r.self_absorption
                    gamma = GetEscapeFactor(r, species_list, system)
                    r.K_value *= gamma
                end

                # Check that rate coeff. is positive
                errcode = K_low_bound_threshold_check(r, system)
                if errcode == c_io_error
                    special_coll_str = ""
                    if special_coll
                        special_coll_str = "special collisions"
                    end
                    message = @sprintf("Interrupt in UpdateRateCoefficientValues! %s", special_coll_str)
                    PrintErrorMessage(system, message)
                    return errcode
                end

            else
                continue
            end
        end

        # Updates regular reactions
        if !special_coll
            if system.prerun
                if r.case == r_extended
                    r.K_value = r.rate_coefficient(dens, temp, species_list,
                        system, sID) 
                else
                    r.K_value = r.rate_coefficient(temp, sID) 
                end
            else
                r.K_value = ReplaceExpressionValues(r.rate_coefficient, temp,
                    species_list, system, sID)
            end
        end

        # Check that rate coeff. is positive
        errcode = K_low_bound_threshold_check(r, system)
        if errcode == c_io_error
            special_coll_str = ""
            if special_coll
                special_coll_str = "special collisions"
            end
            message = @sprintf("Interrupt in UpdateRateCoefficientValues! %s", special_coll_str)
            PrintErrorMessage(system, message)
            return errcode
        end
    end
    return errcode
end


function K_low_bound_threshold_check(r::Reaction, system::System)

    errcode = 0
    if r.K_value < 0.0
        if r.case == r_lower_threshold
            r.K_value = 0.0
        else
            err_message = @sprintf("%i: %s has negative rate coefficient",r.id, r.name)
            PrintErrorMessage(system, err_message)
            return c_io_error
        end
    end
    return errcode

end


function UpdateTotalPressure!(system::System, species_list::Vector{Species}, sID::SpeciesID)
    p_total = 0.0
    for s in species_list
        if s.id == sID.electron
            continue
        end
        p_total += s.pressure 
    end
    system.total_pressure = p_total
    return 0
end


function GetMFP!(species::Species, dens::Vector{Float64}, species_list::Vector{Species})
    # The mean-free-path is calculated using the reaction list associated to
    # species. This reaction_list is attached when reading the input deck and
    # for the case of neutrals and pos/neg-ions only reactions are included which
    # do not involve electrons, i.e. neutral-neutral or neutral-ion reactions
    errcode = 0

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
        cross_section = r.K_value / v_th

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

    species.mfp = 1.0/ilambda
    return errcode 
end


function UpdateElectronegativity!(system::System, dens::Vector{Float64},
    species_list::Vector{Species}, electron_id::Int64)
    # Alpha parameter: is the ratio of negative ion species to electron density
    errcode = 0

    negative_ion_dens = 0.0
    for s in species_list
        if s.id == electron_id
            continue
        elseif s.charge < 0.0
            n = dens[s.id]
            if n > 0.0
                negative_ion_dens += n 
            end
        end
    end

    system.alpha = negative_ion_dens / dens[electron_id]
    return errcode 
end


function GetBohmSpeed!(species::Species, Te::Float64)

    errcode = 0
    try
        species.v_Bohm = sqrt(kb * Te / species.mass)
    catch
        errcode = c_io_error
    end
    return errcode 
end


function GetThermalSpeed!(species::Species)
    errcode = 0

    T = species.temp
    m = species.mass
    try
        species.v_thermal = sqrt(8.0* kb * T /  m / pi)
    catch
        errcode = c_io_error
    end
    return errcode 
end


function GetLambda(system::System)

    L = system.l
    R = system.radius
    Lambda = 1.0/sqrt((pi/L)^2 + (2.405/R)^2)

    return Lambda
end


function GetNeutralDiffusionCoeff!(species::Species)
    # Neutral Diffusion Coefficient

    mfp = species.mfp
    species.D = e * species.temp * mfp / (species.v_thermal * species.mass)

    return 0 
end


function DiffusionRateCoefficient(species::Species, system::System)
    V = system.V
    A = system.A
    lambda = system.Lambda
    D = species.D
    gamma = species.gamma
    vth = species.v_thermal

    K = 1.0/(lambda^2 / D  + 2.0 * V * (2.0 - gamma) / A / vth  / gamma)
    return K
end


function GetSheathDensity!(species::Species, species_list::Vector{Species},
    system::System, sID::SpeciesID)
    # Calculates the density at sheath edge
    # This is only interesting for positive ions. Used later on to
    # calculate the wall-flux = u_B * n_Sheath
    errcode = 0

    n_0 = species.dens

    if system.h_id == h_classical
        # Electropositive plasmas
        L = system.l
        mfp = species.mfp

        if system.total_pressure < 0.26664474 # ( = 2 mTorr) 
            # Low pressure range
            h = 0.425 
        elseif system.total_pressure < 20 # ( = 150 mTorr)
            # Intermediate pressure range
            h = 0.86 / sqrt(3.0 + L/mfp)
        else
           # High pressure
           uB = species.v_Bohm 
           uTh = species.v_thermal 
           h = pi * (uB / uTh) * (mfp / L)
        end
    elseif system.h_id == h_Gudmundsson
        # ICP Ar/O2
        R = system.radius
        L = system.l
        Te = species_list[sID.electron].temp
        Ti = species.temp
        gamma = Te / Ti
        alpha = system.alpha
        lambda = species.mfp
        D = species.D
        u_B = species.v_Bohm

        low_press_term  = 3.0
        int_press_term = 0.5*L/lambda
        high_press_term = 0.86 * L * u_B / (pi * D)
        fL = 1.0/ sqrt(low_press_term + int_press_term + high_press_term^2)
        h_L = 0.86 * (1.0 + 3.0 * alpha/gamma)/(1+gamma) * fL

        low_press_term  = 4.0
        int_press_term = R/lambda
        chi01 = 2.405
        J1_chi01 = 0.52
        high_press_term = 0.8 * R * u_B / (chi01 * J1_chi01 * D)
        fR = 1.0 / sqrt(low_press_term + int_press_term + high_press_term^2)
        h_R = 0.80 * (1.0 + 3.0 * alpha/gamma)/(1+gamma) * fR

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
            K_recombination = GetRecombinationRate(species)
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
            errcode = c_io_error
            PrintErrorMessage(system, "No sheath model found") 
        end
    end
    species.n_sheath = n_0 * h

    return errcode 
end


function GetRecombinationRate(species::Species) 

    s_id = species.id

    K_recomb= 0.0
    for r in species.reaction_list
        s_index = findall( x -> x == s_id, r.involved_species )[1]
        sign = r.species_balance[s_index]
        if sign < 0
            K_recomb += r.K_value
        end
    end
    return K_recomb
end


function GetStickingCoefficient!(species::Species,
    species_list::Vector{Species}, sID::SpeciesID)

    errcode = 0

    if species.id == sID.O
        # This is for stainless steel walls ( Gudmundsson 2007)
        s = species_list[sID.O2]
        pO2_mTorr = s.dens * kb * s.temp * 0.13332237

        if pO2_mTorr < 2.0
            gamma = 1.0 - pO2_mTorr * 0.25
        else
            gamma = 0.1438 * exp(2.5069/pO2_mTorr)
        end
        s.gamma = gamma
    end

    return errcode 
end


function GetEscapeFactor(reaction::Reaction, species_list::Vector{Species},
    system::System)

    # Selected species: absorbing species (lower energy state)
    species = species_list[reaction.product_species[1]]

    # Absorption coefficient
    kappa = GetAbsorptionCoefficient(reaction, species)
    L = 0.5 * ( system.l + 2.0 * system.radius)
    kappa_L = kappa * L 

    # Escape factor
    gamma = (2.0 - exp(- kappa_L * 1.e-3)) / (1.0 + kappa_L)
    gamma *= reaction.g_high /reaction.g_high_total

    return gamma
end


function GetAbsorptionCoefficient(reaction::Reaction, species::Species)

    wavelen2 = reaction.wavelength * reaction.wavelength
    g_rate = reaction.g_high / reaction.g_low
    P_pk = GetSpectralLineProfile_DopplerBroadening(reaction, species, 0.0)
    A_pk = reaction.K_value 
    dens = maximum([0.0,species.dens]) * reaction.g_low / reaction.g_low_total
    kappa = wavelen2 / (8.0 * pi) * P_pk * g_rate * dens * A_pk 
    return kappa
end


function GetSpectralLineProfile_DopplerBroadening(reaction::Reaction,
    species::Species, freq_broadening::Float64)

    wavelen = reaction.wavelength
    wavelen2 = wavelen * wavelen
    mass = species.mass
    temp = species.temp
    ivel2_fac = mass/(2.0*kb*temp)

    P = wavelen * sqrt( ivel2_fac/pi ) *
        exp(-wavelen2 * ivel2_fac * freq_broadening*freq_broadening)
    return P
end

end
