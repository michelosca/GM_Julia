module PlasmaParameters

using SharedData: Species, Reaction, System
using SharedData: kb
using InputBlock_Species: s_electron_id
using InputBlock_Reactions: r_wall_loss

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - GetBohmSpeed
# - GetThermalSpeed
# - GetMFP

function GetMFP(dens::Vector{Float64}, temp::Vector{Float64},
    species::Species, species_list::Vector{Species},
    reaction_list::Vector{Reaction})
    # Get the species-neutral mean-free-path
    # Input species must have non-neutral electric charge

    ilambda = 0.0         # inverse mean-free-path
    s_id = species.id     # species id
    v_th_s = GetThermalSpeed(temp, species)
    for s in species_list
        # Iterate over all species. If neutral species -> check whether there
        # are collisions between s and n species
        if !(s.charge == 0)
            #exclude charged species
            continue
        else
            n_id = s.id    # neutral species id
            n_n = dens[n_id]
            sigma_sn = 0.0 # ion-neutral collision cross-section

            # Loop over reaction set and find reactions where i and n species are involved
            for r in reaction_list
                if r.id == r_wall_loss
                    continue
                end
                # is species s involved?
                i_involved = findall(x->x==s_id, r.reactant_species)
                if i_involved==Int64[]
                    continue
                end

                # is species n involved?
                n_involved = findall(x->x==n_id, r.reactant_species)
                if n_involved==Int64[]
                    continue
                end

                # if here, s and n species are involved -> work out cross-section value
                sigma_sn += r.rate_coefficient(temp) / v_th_s
            end

            ilambda += n_n * sigma_sn
        end
    end
    if ilambda == 0
        print("MFP of species ",s_id," is infinite!\n")
        ilambda = 1.e-100
    end
    lambda = 1.0/ilambda
    return lambda
end


function GetBohmSpeed(temp::Vector{Float64}, species::Species)
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

end