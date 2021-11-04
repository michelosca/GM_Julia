module WallFlux

using SharedData: System, Species, Reaction
using SharedData: kb, e, eps0
using SharedData: p_ccp_id, p_icp_id 
using InputBlock_Species: s_electron_id, s_OnIon_id
using InputBlock_Reactions: r_wall_loss
using Roots: find_zeros

###############################################################################
# WALL FLUX EXPRESSIONS
#   - GetBohmVelocity
#   - GetThermalSpeed
#   - GetIonNeutralMFP
#   - WallFluxFunction_ICP
#   - InterpolateSheathVoltage
#   - MeanSheathVoltage

function GetParticleFlux(dens::Vector{Float64}, temp::Vector{Float64},
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System, V_sheath::Float64)

    flux_list = zeros(length(species_list))
    if (system.power_input_method == p_ccp_id)
        flux_electrons = 0.0
    end
    for s in species_list
        if s.charge == 0
            continue
        end

        if (system.power_input_method == p_ccp_id)
            if s.id == s_electron_id
                continue
            else
                flux_s = ParticleFlux_CCP(dens, temp, s, species_list,
                    reaction_list, system)
                flux_electrons = flux_s * s.charge
            end
        elseif (system.power_input_method == p_icp_id)
            flux_s = ParticleFlux_ICP(dens, temp, s, species_list,
                reaction_list, system, V_sheath)
        end
        flux_list[s.id] = flux_s
    end

    if (system.power_input_method == p_ccp_id)
        e_species = species_list[s_electron_id]
        flux_list[e_species.id] = flux_electrons / abs(e_species.charge)
    end

    return flux_list
end


function ParticleFlux_CCP(dens::Vector{Float64},
    temp::Vector{Float64},
    species::Species,
    species_list::Vector{Species},
    reaction_list::Vector{Reaction},
    system::System)

    ni = dens[species.id] # Ion density
    uB = GetBohmVelocity(temp, species) # Bohm velocity
    uTh = GetThermalSpeed(temp, species) # Thermal speed
    lambda = GetIonNeutralMFP(dens, temp, species, species_list, reaction_list) # MFP
    l = system.l # System length
    #print("  - Species ", species.id, " ni:",ni," uB:", uB," uTh:", uTh," mfp:",lambda,"\n")

    n_sheath = pi * ni * uB / uTh * lambda / l # Density at sheath edge
    wallflux = n_sheath * uB
    #print("  - Species ", species.id, " n_sheath:",n_sheath," flux:", wallflux,"\n")

    return wallflux
end


function ParticleFlux_ICP(dens, temp, species, species_list, reaction_list,
    system, V_sheath)

    charge = species.charge
    if charge > 0
        R = system.radius
        L = system.l
        Ti = temp[species.id]
        Te = temp[s_electron_id]
        if s_OnIon_id == 0
            alpha = 0.0
        else
            ne = dens[s_electron_id]
            n0neg = dens[s_OnIon_id]
            alpha = n0neg / ne
        end
        gamma = Te / Ti
        lambda_i = GetIonNeutralMFP(dens, temp, species, species_list, reaction_list)
        h_L = 0.86 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(3.0+0.5*L/lambda_i)
        h_R = 0.8 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(4.0+R/lambda_i)
        fac = (R^2 * h_L + R*L*h_R) / (R^2 + R*L)
        uB_i = GetBohmVelocity(temp, species)
        n_i = dens[species.id]
        flux = uB_i * fac * n_i
        return flux
    elseif charge < 0
        T_eV = temp[species.id]/K_to_eV
        n = dens[species.id]
        vth = GetThermalSpeed(temp, species)
        flux_funct = (pot) -> 0.25 * n * vth *exp(-pot/T_eV)
        if V_sheath == 0
            # return potential functions
            return flux_funct
        else
            return flux_funct(V_sheath)
        end
    end
    return 0
end


function GetIonNeutralMFP(dens::Vector{Float64}, temp::Vector{Float64},
    ion_species::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction})
    # Get the ion-neutral mean-free-path
    # Input species must be a positive ion

    if (ion_species.charge <= 0)
        print("Lambda values only make sense with positive ions!\n")
        return 0
    end

    ilambda = 0.0         # inverse mean-free-path
    i_id = ion_species.id # ion lambda species
    v_th_i = GetThermalSpeed(temp, ion_species)
    for s in species_list
        # Iterate over all species. If neutral species -> check whether there
        # are collisions between i and j species
        if !(s.charge == 0)
            #exclude charged species
            continue
        else
            n_id = s.id    # neutral species id
            n_n = dens[n_id]
            sigma_in = 0.0 # ion-neutral collision cross-section

            # Loop over reaction set and find reactions where i and n species are involved
            for r in reaction_list
                if r.id == r_wall_loss
                    continue
                end
                # is species i involved?
                i_involved = findall(x->x==i_id, r.reactant_species)
                if i_involved==Int64[]
                    continue
                end

                # is species n involved?
                n_involved = findall(x->x==n_id, r.reactant_species)
                if n_involved==Int64[]
                    continue
                end

                # if here, i and n species are involved -> work out cross-section value
                sigma_in += r.rate_coefficient(temp) / v_th_i
                #print("  - ** Species ", i_id, "-", n_id, " reaction:",r.id," sigma:", sigma_in,"\n")
            end

            ilambda += n_n * sigma_in
        end
    end
    if ilambda == 0
        print("MFP of species ",i_id," is infinite!\n")
        ilambda = 1.e-100
    end
    lambda = 1.0/ilambda
    return lambda
end

function GetBohmVelocity(temp::Vector{Float64}, species::Species)
    if (species.id == s_electron_id)
        print("Electron Bohm velocity does not make sense!\n")
        return 0
    end
    Te = temp[s_electron_id]
    ms = species.mass
    uB = sqrt(kb * Te / ms)
    return uB
end

function GetThermalSpeed(temp::Vector{Float64}, species::Species)

    T = temp[species.id]
    m = species.mass
    v_th = sqrt(8.0* kb * T /  m / pi)
    return v_th
end

###############################################################################
### SHEATH FUNCTIONS

function MeanSheathVoltage(dens::Vector{Float64}, system::System)

    ne = dens[s_electron_id]
    V_sheath = 3.0/4.0 * (system.drivI/system.A)^2 / (e * eps0 * ne * system.drivOmega^2)
    return V_sheath
end


function InterpolateSheathVoltage(dens, temp, species_list,
    reaction_list, system)
    # Obtain the positive and negatice charged particles flux and solve for the potential

    flux_plus = 0.0
    flux_min = Function[]
    for s in species_list
        if s.charge > 0
            flux_plus += ParticleFlux_ICP(dens, temp, s, species_list,
                reaction_list, system, 0.0)
        elseif s.charge < 0
            push!(flux_min, ParticleFlux_ICP(dens, temp, s, species_list,
                reaction_list, system, 0.0))
        end
    end

    # Solve now "flux_min - flux_plus = 0" and get the plasma potential
    function flux_min_funct(pot)
        output = 0.0
        for f in flux_min
            output += f(pot)
        end
        return output - flux_plus
    end
    V_sheath = find_zeros(flux_min_funct, 0)

    return V_sheath
end

###############################################################################
### Density and temperature rate fluxes

function DensWallFluxFunction(flux_list::Vector{Float64}, system::System,
    species::Species)
    A = system.A
    V = system.V
    dens_flux = -A/V * flux_list[species.id]
    return dens_flux
end

function TempWallFluxFunction(temp::Vector{Float64}, flux_list::Vector{Float64},
    system::System, V_sheath::Float64, species::Species)
    A = system.A
    V = system.V
    Te = temp[s_electron_id]
    if species.id == s_electron_id
        fact = 2.0 * kb * Te
    else
        fact = 0.5 * kb * Te + abs(s.charge) * V_sheath
    end
    temp_flux = -fact*A/V * flux_list[species.id]
    return temp_flux
end

end