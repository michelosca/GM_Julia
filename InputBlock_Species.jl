module InputBlock_Species

using SharedData: K_to_eV, e, me, amu 
using SharedData: c_io_error
using SharedData: Species, Reaction, SpeciesID

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
    if (read_step == 2)
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
        current_species.reaction_list = Reaction[]
        current_species.mfp = 1.e100
        current_species.v_thermal = 0.0
        current_species.v_Bohm = 0.0
        current_species.D = 0.0
        current_species.Lambda = 0.0
        current_species.h_R = 0.0
        current_species.h_L = 0.0
        current_species.gamma = 0.0
        current_species.n_sheath = 0.0
        current_species.flux = 0.0
        push!(species_list, current_species)
    end
    return errcode
end


function ReadSpeciesEntry!(name::SubString{String}, var::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, sID::SpeciesID)

    errcode = c_io_error
    if (read_step == 2)
        current_species = species_list[end]
    end

    if (name=="name")
        if (read_step == 2)
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
        else 
            # This is set in the pre-run(read_step==0) and in main-run(read_step===1)
            errcode = SetSpeciesID!(var, sID)
        end
    end

    if (name=="charge")
        if (read_step == 2)
            current_species.charge = parse(Int64,var) * e
        end
        errcode = 0
    end

    if (name=="mass")
        if (read_step == 2)
            expr = Meta.parse(var)
            current_species.mass = eval(expr)
        end
        errcode = 0
    end

    if (name=="solve_dens")
        if (read_step == 2)
            current_species.has_dens_eq = parse(Bool, var) 
        end
        errcode = 0
    end

    if (name=="solve_temp")
        if (read_step == 2)
            current_species.has_temp_eq = parse(Bool, var)
        end
        errcode = 0
    end

    if (name=="T" || name=="T_eV")
        if (read_step == 2)
            if (name=="T_eV")
                units = 1.0/K_to_eV
            else
                units = 1.0
            end
            current_species.temp = parse(Float64, var)*units
        end
        errcode = 0
    end

    if (name=="density" || name=="dens")
        if (read_step == 2)
            current_species.dens = parse(Float64, var)
        end
        errcode = 0
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
    elseif ("O(3p)" == species_name)
        speciesID.O_3p = id
    elseif ("O(1d)" == species_name)
        speciesID.O_1d = id
    elseif ("O2(a1Ag)" == species_name)
        speciesID.O2_a1Ag = id
    else
        errcode = 1
    end
    return errcode 
end


function EndSpeciesBlock!(read_step::Int64, species_list::Vector{Species})

    errcode = 0 

    if (read_step == 1)
        errcode = 0
    elseif (read_step == 2)
        current_species = species_list[end]
        if (current_species.charge != 0)
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
            print("***ERROR*** Species temperature has not been defined")
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
    speciesID.O_3p = 0

    speciesID.O2 = 0
    speciesID.O2_Ion = 0
    speciesID.O2_a1Ag = 0

end

end