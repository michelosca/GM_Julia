module WallFlux

using SharedData: System, Species, Reaction
using SharedData: kb, e 
using SharedData: p_ccp_id, p_icp_id 
using PlasmaParameters: GetMFP, GetBohmSpeed, GetThermalSpeed
using InputBlock_Species: s_electron_id, s_OnIon_id

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - DensWallFluxFunction
# - TempWallFluxFunction
# - GetIonFlux
#   - IonParticleFlux_ICP
#   - IonParticleFlux_CCP
# - GetElectronFlux!

function DensWallFluxFunction(flux_list::Vector{Tuple}, system::System,
    species::Species)
    A = system.A
    V = system.V
    dens_flux = -A/V * flux_list[species.id][1]
    return dens_flux
end

function TempWallFluxFunction(temp::Vector{Float64}, flux_list::Vector{Tuple},
    species_list::Vector{Species} ,system::System, V_sheath::Float64, species::Species)
    # input "species" is the species whose temperature equation this term is going to be
    A = system.A
    V = system.V
    s_id = species.id
    if s_id == s_electron_id
        Te = temp[s_electron_id]
        fact = 2.0 * kb * Te * flux_list[s_id][1]
        # Add ion flux terms
        for s in species_list
            charge = s.charge
            if (charge > 0)
                fact += (0.5*kb*Te + charge*V_sheath) * flux_list[s.id][1]
            end
        end
    else
        fact = 0.0
        print("Temperature flux model need to be defined for species ", s_id,"\n")
    end
    temp_flux = -fact*A/V
    return temp_flux
end


function GetIonFlux(dens::Vector{Float64}, temp::Vector{Float64},
    species_list::Vector{Species}, reaction_list::Vector{Reaction},
    system::System)

    flux_list = Tuple[]
    for s in species_list
        if s.charge <= 0
            push!(flux_list, (0.0,0.0))
        else
            if (system.power_input_method == p_ccp_id)
                flux_tuple = IonParticleFlux_CCP(dens, temp, s, species_list,
                    reaction_list, system)
                push!(flux_list, flux_tuple)
            elseif (system.power_input_method == p_icp_id)
                flux_s = IonParticleFlux_ICP(dens, temp, s, species_list,
                    reaction_list, system)
                push!(flux_list, (flux_s, 0.0))
            end
        end
    end

    return flux_list
end


function GetElectronFlux!(flux_list::Vector{Tuple},
    species_list::Vector{Species})

    flux_electrons = 0.0
    n_sheath_e = 0.0
    for s in species_list
        if s.charge > 0
            flux_electrons += flux_list[s.id][1] * s.charge
            n_sheath_e += flux_list[s.id][2] * s.charge
        end
    end
    flux_list[s_electron_id] = ( flux_electrons/e, n_sheath_e/e )

    return flux_list
end


function IonParticleFlux_CCP(dens::Vector{Float64},
    temp::Vector{Float64},
    species::Species,
    species_list::Vector{Species},
    reaction_list::Vector{Reaction},
    system::System)

    charge = species.charge
    if charge > 0
        ni = dens[species.id] # Ion density
        uB = GetBohmSpeed(temp, species) # Bohm velocity
        uTh = GetThermalSpeed(temp, species) # Thermal speed
        lambda = GetMFP(dens, temp, species, species_list, reaction_list) # MFP
        l = system.l # System length

        n_sheath = pi * ni * (uB / uTh) * (lambda / l) # Density at sheath edge
        flux = n_sheath * uB
    else
        n_sheath = 0.0
        flux = 0.0
    end
    return (flux, n_sheath)
end


function IonParticleFlux_ICP(dens::Vector{Float64}, temp::Vector{Float64},
    species::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

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
        lambda_i = GetMFP(dens, temp, species, species_list, reaction_list)
        h_L = 0.86 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(3.0+0.5*L/lambda_i)
        h_R = 0.8 * (1.0 + 3.0 * alpha/gamma)/(1+gamma)/sqrt(4.0+R/lambda_i)
        fac = (R^2 * h_L + R*L*h_R) / (R^2 + R*L)
        uB_i = GetBohmSpeed(temp, species)
        n_i = dens[species.id]
        flux = uB_i * fac * n_i
    else
        flux  = 0.0
    end
    return flux
end

end