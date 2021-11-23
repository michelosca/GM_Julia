module PlasmaParameters

using SharedData: Species, Reaction, System, SpeciesID
using SharedData: kb, K_to_eV, e
using SharedData: p_icp_id, p_ccp_id
using InputBlock_Reactions: r_wall_loss
using EvaluateExpressions: EvaluateExpression

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - UpdateSpeciesParameters
# - GetBohmSpeed
# - GetThermalSpeed
# - GetMFP

function UpdateSpeciesParameters!(temp::Vector{Float64}, dens::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)

    for s in species_list
        s.temp = temp[s.id]
        s.dens = dens[s.id]
        s.v_thermal = GetThermalSpeed(s)
        s.v_Bohm = GetBohmSpeed(temp[sID.electron], s.mass)
        s.mfp = GetMFP(temp, dens, s, sID)
        if system.power_input_method == p_ccp_id
            if s.charge > 0
                s.n_sheath = GetSheathDensity(s, system)
            else
                s.n_sheath = s.dens
            end
        elseif system.power_input_method == p_icp_id
            s.n_sheath = s.dens
            if (s.id == sID.electron)
                continue
            end
            #if s.charge > 0
            #    s.n_sheath = GetSheathDensity(s, system)
            #end
            s.h_L, s.h_R = Get_h_Parameters(temp, dens, s, species_list,
                system, sID.electron)
            s.D = GetD(s)
        end
    end
end

function GetMFP(temp::Vector{Float64}, dens::Vector{Float64}, species::Species,
    sID::SpeciesID)
    # Get the species neutral mean-free-path

    ilambda = 0.0         # inverse mean-free-path
    v_th_s = species.v_thermal
    for r in species.reaction_list
        if (r.case == r_wall_loss)
            continue
        end
        # Set collision cross section
        sigma_expr = EvaluateExpression(r.rate_coefficient, sID, temp, species)
        cross_section = eval(sigma_expr) / v_th_s

        # Density of colliding partners
        n = prod(dens[r.neutral_species_id])
        ilambda += n * cross_section 
    end

    if ilambda == 0
        ilambda = 1.e-100
    end
    lambda = 1.0/ilambda
    return lambda
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

    L = system.l
    R = system.radius

    # Gamma parameter
    Te = temp[electron_id]
    Ti = temp[species.id]
    gamma = Te / Ti

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
    D = e * T_eV * species.mfp / (species.v_thermal / species.mass)
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