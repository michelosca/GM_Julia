module WallFluxModels

using SharedData
using SharedData: A, V, l, drivOmega, drivI
using SharedData: kb, e, eps0
using SharedData: s_electron_id
using SharedData: Species

function WallFluxFunction_CCP(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species)
    # Input: ions species
    #        if species is electron then summ of all ion species

    if (species.id == s_electron_id)
        # Electron wall flux = sum of all ion wall fluxes
        wallflux = 0
        for s in SharedData.species_list
            if s.charge>0
                wallflux += WallFluxFunction_CCP_ion_species(dens, temp, s)
            end
        end
    else
        wallflux = WallFluxFunction_CCP_ion_species(dens, temp, species)
    end
    return wallflux
end


function WallFluxFunction_CCP_ion_species(dens::Vector{Float64},
    temp::Vector{Float64}, s::Species)
    # Input: Species type (must be an ion!)
    
    # Plasma parameters 
    Te = temp[s_electron_id]
    ni = dens[s.id]
    mi = s.mass

    # Ion total collision frequency (elastic collision)
    r_elastic = SharedData.reaction_list[s.r_elastic_id]
    ng = dens[r_elastic.neutral_species_id]
    K_elastic = r_elastic.rate_coefficient(temp)
    nu = ng * K_elastic

    wallflux = ni * π * kb * Te / l / mi / nu
    return wallflux

end


function MeanSheathVoltage(dens::Vector{Float64})

    ne = dens[s_electron_id]
    V_sheath = 3/4 * (drivI/A)^2 / (e * eps0 * ne * drivOmega^2)
    return V_sheath
end


function TempWallFluxFunction(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species)
    # Input: Species type (must be an ion!)

    if (species.charge == 0)
        temp_wf = 0
    else
        if (SharedData.power_input_method == SharedData.p_ccp_id)
            if (species.id == s_electron_id)
                # Must account for the electron flux and for the ion flux
                temp_wf = 0
                V_sheath = MeanSheathVoltage(dens)
                Te = temp[s_electron_id]
                for s in SharedData.species_list
                    if (s.id == s_electron_id)
                        # Electron flux
                        wf = WallFluxFunction_CCP(dens, temp, s)
                        Te = temp[s_electron_id]
                        temp_wf -= wf * A / V * 2 * kb * Te
                    elseif (s.charge>0)
                        # Ion flux
                        q = s.charge
                        wf = WallFluxFunction_CCP(dens, temp, s)
                        temp_wf -= wf * A / V * (0.5*kb*Te + q*V_sheath)
                    end
                end
            else
                print("Temperature flux models for ions have not been implemented\n")
                temp_wf = 0
            end
        else
            print("Flux models only implemented for CCPs\n")
            temp_wf = 0
        end
    end
    return temp_wf

end

function DensWallFluxFunction(dens::Vector{Float64},
    temp::Vector{Float64}, species::Species)

    if (species.charge == 0)
        dens_wf = 0
    else
        if (SharedData.power_input_method == SharedData.p_ccp_id)
            wf = WallFluxFunction_CCP(dens, temp, species)
            dens_wf = -A / V * wf
        else
            dens_wf = 0
            print("Flux models only implemented for CCPs\n")
        end
    end
    return dens_wf
end

end