module WallFlux

using SharedData: System, Species, Reaction
using SharedData: kb, e 
using SharedData: p_ccp_id, p_icp_id 
using PlasmaParameters: GetMFP, GetBohmSpeed, GetThermalSpeed
using PlasmaParameters: Get_h_Parameters
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

function DensWallFluxFunction(species::Species)
    return -species.flux 
end

function TempWallFluxFunction(temp::Vector{Float64}, species::Species,
    species_list::Vector{Species}, V_sheath::Float64)

    s_id = species.id
    temp_flux = 0.0

    if s_id == s_electron_id
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

    return temp_flux
end


function UpdateIonFlux!(species_list::Vector{Species}, system::System)

    # Geometrical factors
    if system.power_input_method == p_ccp_id
        V = system.V
        A = system.A
        fac = A/V
    elseif system.power_input_method == p_icp_id
        L = system.l
        R = system.radius
    end

    # Particle flux of positive ions
    for s in species_list
        if s.charge > 0
            if system.power_input_method == p_ccp_id
                s.flux = fac * s.v_Bohm * s.n_sheath 
            elseif system.power_input_method == p_icp_id
                h_L = s.h_L
                h_R = s.h_R
                fac = (R^2 * h_L + R*L*h_R) / (R^2 + R*L)
                s.flux = fac * s.v_Bohm * s.dens
            end
        end
    end
end


function UpdateElectronFlux!(species_list::Vector{Species})

    # Electron flux solved from the flux balance equation, i.e.
    # flux_electrons = SUM(flux_ions)
    # Negative ion species are assumed to have zero net flow through the walls

    flux_electrons = 0.0
    for s in species_list
        if s.charge > 0
            flux_electrons += s.flux * s.charge
        end
    end
    species_list[s_electron_id].flux = flux_electrons / e

end

end