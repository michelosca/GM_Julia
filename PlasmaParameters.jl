module PlasmaParameters

using SharedData: Species, Reaction, System
using SharedData: kb
using SharedData: p_icp_id, p_ccp_id
using InputBlock_Species: s_electron_id
using InputBlock_Reactions: r_wall_loss

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - UpdateSpeciesParameters
# - GetBohmSpeed
# - GetThermalSpeed
# - GetMFP

function UpdateSpeciesParameters!(temp::Vector{Float64}, dens::Vector{Float64},
    species_list::Vector{Species}, system::System)

    for s in species_list
        s.temp = temp[s.id]
        s.dens = dens[s.id]
        s.v_thermal = GetThermalSpeed(s)
        s.v_Bohm = GetBohmSpeed(temp, s)
        s.mfp = GetMFP(temp, dens, s)
        if system.power_input_method == p_ccp_id
            if s.charge > 0
                s.n_sheath = GetSheathDensity(s, system)
            end
        elseif system.power_input_method == p_icp_id
            if (s.id == s_electron_id)
                continue
            end
            if s.charge > 0
                s.n_sheath = s.dens
            end
            s.h_L, s.h_R = Get_h_Parameters(temp, dens, s, species_list, system)
            s.D = GetD(s)
        end
    end
end


function GetMFP(temp::Vector{Float64}, dens::Vector{Float64}, species::Species)
    # Get the species neutral mean-free-path

    ilambda = 0.0         # inverse mean-free-path
    v_th_s = GetThermalSpeed(species)
    for r in species.reaction_list
        cross_section = r.rate_coefficient(temp) / v_th_s
        n = prod(dens[r.reactant_species]) / species.dens
        ilambda += n * cross_section 
    end

    if ilambda == 0
        #print("MFP of species ",species.id," is infinite!\n")
        ilambda = 1.e-100
    end
    lambda = 1.0/ilambda
    return lambda
end


function GetBohmSpeed(temp::Vector{Float64}, species::Species)
    if (species.id == s_electron_id)
        return 0.0
    end
    Te = temp[s_electron_id]
    ms = species.mass
    uB = sqrt(kb * Te / ms)
    return uB
end


function GetThermalSpeed(species::Species)

    T = species.temp
    m = species.mass
    v_th = sqrt(8.0* kb * T /  m / pi)
    return v_th
end


function Get_h_Parameters(temp::Vector{Float64}, dens::Vector{Float64}, species::Species,
    species_list::Vector{Species}, system::System)
    # See Gudmundsson (2000) On the plasma parameters of a planar inductive oxygen discharge

    L = system.l
    R = system.radius

    # Gamma parameter
    Te = temp[s_electron_id]
    Ti = temp[species.id]
    gamma = Te / Ti

    # Alpha parameter
    alpha = 0.0
    for s in species_list
        q = s.charge
        if s.id == s_electron_id
            continue
        elseif q < 0
            alpha += dens[s.id]
        end
    end
    alpha /= dens[s_electron_id]

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


function GetGamma()
    return 1.0
end


function GetD(species)

    T_eV = species.temp * K_to_eV
    D = e * T_eV * species.mfp / (species.v_th / species.mass)
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