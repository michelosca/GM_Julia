module InputBlock_Species

using SharedData: K_to_eV, e, me, amu, kb 
using SharedData: c_io_error, p_icp_id
using SharedData: r_wall_loss
using SharedData: Species, Reaction, SpeciesID, System
using PlasmaParameters: GetGamma
using InputBlock_System: GetUnits!

###############################################################################
################################  VARIABLES  ##################################
###############################################################################

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - StartFile_Species
# - StartSpeciesBlock
# - ReadSpeciesEntry
#   - SetSpeciesID
# - EndSpeciesBlock
#   - AddSpeciesToList


function StartFile_Species!(read_step::Int64, species_list::Vector{Species},
    speciesID::SpeciesID) 

    errcode = 0

    if (read_step == 1)
        InitializeSpeciesID!(speciesID)
    end
    speciesID.current_id = 0

    return errcode
end

function StartSpeciesBlock!(read_step::Int64, species_list::Vector{Species},
    speciesID::SpeciesID)

    errcode = 0
    speciesID.current_id += 1
    if (read_step == 1)
        current_species = Species()
        current_species.id = speciesID.current_id 
        current_species.species_id = 0
        current_species.mass = 0.0
        current_species.charge = 0.0
        current_species.has_dens_eq = false
        current_species.has_temp_eq = false
        current_species.has_wall_loss = false
        current_species.has_heating_mechanism = false
        current_species.dens = 0.0
        current_species.temp = 0.0
        current_species.pressure = 0.0
        current_species.reaction_list = Reaction[]
        current_species.mfp = 1.e100
        current_species.cross_section = 0.0
        current_species.v_thermal = 0.0
        current_species.v_Bohm = 0.0
        current_species.D = 0.0
        current_species.h_R = 0.0
        current_species.h_L = 0.0
        current_species.gamma = 0.0
        current_species.n_sheath = 0.0
        current_species.flux = 0.0
        current_species.name = "None"
        current_species.has_flow_rate = false
        current_species.flow_rate = 0.0 
        push!(species_list, current_species)
    end
    return errcode
end


function ReadSpeciesEntry!(name::SubString{String}, var::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, sID::SpeciesID)

    errcode = 0 

    if (read_step == 1)
        units, name = GetUnits!(name)

        if (name=="name")
            # This is set in the pre-run(read_step==0) and in main-run(read_step===1)
            errcode = SetSpeciesID!(var, sID)
            current_species = species_list[sID.current_id]
            current_species.name = var
        end

        current_species = species_list[sID.current_id]

        if (name=="charge")
            current_species.charge = parse(Int64,var) * e * units
        end

        if (name=="mass")
            expr = Meta.parse(var)
            current_species.mass = eval(expr) * units
        end

        if (name=="solve_dens")
            current_species.has_dens_eq = parse(Bool, var) 
        end

        if (name=="solve_temp")
                current_species.has_temp_eq = parse(Bool, var)
        end

        if (name=="T")
            current_species.temp = parse(Float64, var) * units
        end

        if (name=="density" || name=="dens")
            current_species.dens = parse(Float64, var) *units
        end

        if (name=="pressure" || name=="p")
            current_species.pressure = parse(Float64, var) * units
        end

        if (name=="cross_section")
            current_species.cross_section = parse(Float64, var) * units
        end

        if (name=="flow_rate")
            current_species.has_flow_rate = true
            current_species.flow_rate = parse(Float64, var) * units
        end

    elseif (read_step == 2)

        current_species = species_list[sID.current_id]

        if (name=="name")
            if (var=="e" || var=="electrons" || var=="electron")
                current_species.species_id = sID.electron 
                errcode = 0
            elseif (occursin("Ar",var))
                current_species.species_id = sID.Ar
                errcode = 0
            elseif (occursin("O2",var))
                current_species.species_id = sID.O2
                errcode = 0
            elseif (occursin("O",var))
                current_species.species_id = sID.O
                errcode = 0
            else
                print("***ERROR*** Neutral species id has not been found\n")
                errcode = c_io_error 
            end
        end
    end

    return errcode 
end


function SetSpeciesID!(species_name::SubString{String}, speciesID::SpeciesID)

    errcode = 0
    id = speciesID.current_id

    if ("e" == species_name || "electrons" == species_name)
        speciesID.electron = id
    elseif ("Ar" == species_name)
        speciesID.Ar = id
    elseif ("O" == species_name)
        speciesID.O = id
    elseif ("O2" == species_name)
        speciesID.O2 = id
    elseif ("Ar+" == species_name)
        speciesID.Ar_Ion = id
    elseif ("O+" == species_name)
        speciesID.O_Ion = id
    elseif ("O-" == species_name)
        speciesID.O_negIon = id
    elseif ("O2+" == species_name)
        speciesID.O2_Ion = id
    elseif ("Ar*" == species_name || "Ar_excited" == species_name)
        speciesID.Ar_Exc = id
    elseif ("O(1d)" == species_name)
        speciesID.O_1d = id
    elseif ("O2(a1Ag)" == species_name)
        speciesID.O2_a1Ag = id
    else
        errcode = 1
    end
    return errcode 
end


function EndSpeciesBlock!(read_step::Int64, species_list::Vector{Species},
    sID::SpeciesID)

    errcode = 0 

    if (read_step == 1)
        # Update equations/ wall losses flags
        current_species = species_list[end]
        if (current_species.charge > 0 || current_species.id == sID.electron)
            if (current_species.has_temp_eq)
                current_species.has_heating_mechanism = true
            end
            current_species.has_wall_loss = true
        end
        
        if (current_species.mass == 0)
            print("***ERROR*** Species mass has not been defined")
            return c_io_error
        end

        if (current_species.temp == 0)
            dens = current_species.dens
            p = current_species.pressure
            if (dens != 0 && p != 0)
                current_species.temp = p / (kb * dens) 
            end
        end

        if (current_species.dens == 0)
            temp = current_species.temp
            p = current_species.pressure
            if (temp != 0 && p != 0)
                current_species.dens = p / (kb * temp) 
            end
        end

        if (current_species.pressure == 0)
            temp = current_species.temp
            dens = current_species.dens
            if (temp != 0 && dens != 0)
                current_species.pressure = dens * kb * temp 
            end
        end

        if (current_species.temp == 0 && current_species.dens == 0)
            print("***ERROR*** Species temp/dens/press has not been defined\n")
            return c_io_error
        end
    end
    return errcode 
end

function InitializeSpeciesID!(speciesID::SpeciesID)

    speciesID.electron = 0
    
    speciesID.Ar = 0
    speciesID.Ar_Ion = 0
    speciesID.Ar_Exc = 0
    
    speciesID.O = 0
    speciesID.O_negIon = 0
    speciesID.O_Ion = 0
    speciesID.O_1d = 0

    speciesID.O2 = 0
    speciesID.O2_Ion = 0
    speciesID.O2_a1Ag = 0

end


function EndFile_Species!(read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System)

    errcode = 0
    
    if (read_step == 2)

        for s in species_list
            s_id = s.id

            # Create reaction list associated to species s
            for r in reaction_list
                if r.case == r_wall_loss
                    continue
                end
                # is species s involved?
                i_involved = findall(x->x==s_id, r.reactant_species)
                if i_involved==Int64[]
                    continue
                else
                    push!(s.reaction_list, r)
                end
            end

            # Sticking coefficient
            if system.power_input_method == p_icp_id
                s.gamma = GetGamma()
            end
        end

    end
    return errcode
end

end