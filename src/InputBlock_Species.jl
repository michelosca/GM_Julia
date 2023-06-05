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

module InputBlock_Species

using SharedData: K_to_eV, e, me, amu, kb 
using SharedData: c_io_error
using SharedData: r_diffusion, r_emission_rate
using SharedData: Species, Reaction, SpeciesID, System
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
        current_species.in_nodes = Tuple{Int64, String, Float64}[]
        current_species.out_nodes = Tuple{Int64, String, Float64}[]
        push!(species_list, current_species)
    end
    return errcode
end


function ReadSpeciesEntry!(name::SubString{String}, var::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, system::System, sID::SpeciesID)

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
            if system.total_pressure > 0
                pres_var = parse(Float64, var)
                if pres_var >= 0 && pres_var <= 1.0
                    current_species.pressure = pres_var * system.total_pressure
                else
                    print("***ERROR*** Partial pressure is out of range\n")
                    errcode = c_io_error
                end
            else
                current_species.pressure = parse(Float64, var) * units
            end
        end

        if (name=="flow_rate")
            current_species.has_flow_rate = true
            current_species.flow_rate = parse(Float64, var) * units
        end

        if (name=="gamma" || name=="sticking_coefficient")
            current_species.gamma = parse(Float64, var)
        end
    end

    return errcode 
end


function SetSpeciesID!(species_name::SubString{String}, speciesID::SpeciesID)

    errcode = 0
    id = speciesID.current_id

    # ELECTRONS
    if ("e" == species_name || "electrons" == species_name)
        speciesID.electron = id

    # ARGON
    elseif ("Ar" == species_name)
        speciesID.Ar = id
    elseif ("Ar+" == species_name)
        speciesID.Ar_Ion = id
    elseif ("Ar_m" == species_name)
        speciesID.Ar_m = id
    elseif ("Ar_r" == species_name)
        speciesID.Ar_r = id
    elseif ("Ar_4p" == species_name)
        speciesID.Ar_4p = id

    # ATOMIC OXYGEN
    elseif ("O" == species_name)
        speciesID.O = id
    elseif ("O+" == species_name)
        speciesID.O_Ion = id
    elseif ("O-" == species_name)
        speciesID.O_negIon = id
    elseif ("O_1d" == species_name)
        speciesID.O_1d = id
    elseif ("O_1s" == species_name)
        speciesID.O_1s = id
    elseif ("O_3s" == species_name)
        speciesID.O_3s = id
    elseif ("O_5s" == species_name)
        speciesID.O_5s = id
    elseif ("O_3p" == species_name)
        speciesID.O_3p = id
    elseif ("O_5p" == species_name)
        speciesID.O_5p = id

    # MOLECULAR OXYGEN
    elseif ("O2" == species_name)
        speciesID.O2 = id
    elseif ("O2_v" == species_name)
        speciesID.O2_v = id
    elseif ("O2+" == species_name)
        speciesID.O2_Ion = id
    elseif ("O2-" == species_name)
        speciesID.O2_negIon = id
    elseif ("O2_a1Ag" == species_name)
        speciesID.O2_a1Ag = id
    elseif ("O2_b1Su" == species_name)
        speciesID.O2_b1Su = id
    elseif ("O2_a1Ag_v" == species_name)
        speciesID.O2_a1Ag_v = id
    elseif ("O2_b1Su_v" == species_name)
        speciesID.O2_b1Su_v = id

    # OZONE and O4
    elseif ("O3" == species_name)
        speciesID.O3 = id
    elseif ("O3_v" == species_name)
        speciesID.O3_v = id
    elseif ("O3+" == species_name)
        speciesID.O3_Ion = id
    elseif ("O3-" == species_name)
        speciesID.O3_negIon = id
    elseif ("O4" == species_name)
        speciesID.O4 = id
    elseif ("O4+" == species_name)
        speciesID.O4_Ion = id
    elseif ("O4-" == species_name)
        speciesID.O4_negIon = id

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
        if !(current_species.charge == 0)
            current_species.has_wall_loss = true
        end

        if current_species.has_temp_eq
            current_species.has_heating_mechanism = true
        end
        
        if current_species.mass == 0
            print("***ERROR*** Species mass has not been defined")
            return c_io_error
        end

        if current_species.temp == 0
            dens = current_species.dens
            p = current_species.pressure
            if (dens != 0 && p != 0)
                current_species.temp = p / (kb * dens) 
            end
        end

        if current_species.dens == 0
            temp = current_species.temp
            p = current_species.pressure
            if (temp != 0 && p != 0)
                current_species.dens = p / (kb * temp) 
            end
        end

        if current_species.pressure == 0
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
    speciesID.Ar_m = 0
    speciesID.Ar_r = 0
    speciesID.Ar_4p = 0
    
    speciesID.O = 0
    speciesID.O_negIon = 0
    speciesID.O_Ion = 0
    speciesID.O_1d = 0
    speciesID.O_1s = 0
    speciesID.O_3s = 0
    speciesID.O_5s = 0
    speciesID.O_3p = 0
    speciesID.O_5p = 0

    speciesID.O2 = 0
    speciesID.O2_v = 0
    speciesID.O2_Ion = 0
    speciesID.O2_negIon = 0
    speciesID.O2_a1Ag = 0
    speciesID.O2_a1Ag_v = 0
    speciesID.O2_b1Su = 0
    speciesID.O2_b1Su_v = 0

    speciesID.O3 = 0
    speciesID.O3_v = 0
    speciesID.O3_Ion = 0
    speciesID.O3_negIon = 0
    speciesID.O4 = 0
    speciesID.O4_Ion = 0
    speciesID.O4_negIon = 0
end


function EndFile_Species!(read_step::Int64, species_list::Vector{Species},
    reaction_list::Vector{Reaction}, system::System, sID::SpeciesID)

    errcode = 0
    
    if (read_step == 2)

        id_electrons = sID.electron

        if system.total_pressure > 0
            pressure_check = 0
        end

        for s in species_list
            s_id = s.id

            # Set the "characteristic" species ID
            var = s.name
            if (var=="e" || var=="electrons" || var=="electron")
                s.species_id = sID.electron 
            elseif (occursin("Ar",var))
                s.species_id = sID.Ar
            elseif (occursin("O4",var))
                s.species_id = sID.O4
            elseif (occursin("O3",var))
                s.species_id = sID.O3
            elseif (occursin("O2",var))
                s.species_id = sID.O2
            elseif (occursin("O",var))
                s.species_id = sID.O
            end
            if s.species_id == 0 
                print("***ERROR*** Neutral species id has not been found\n")
                errcode = c_io_error 
            end

            # Create reaction list associated to species s
            # This is ONLY used in calculating the MFP, therefore
            # - for ION SPECIES (positive and negative) only ion-neutral
            #  reactions are included
            # - for NEUTRAL SPECIES only ion-neutral reactions are included
            # - for ELECTRONS any reaction is included
            for r in reaction_list
                if r.case == r_diffusion || r.case == r_emission_rate
                    continue
                end

                # is species s involved?
                if (s_id == id_electrons)
                    # For the electron species
                    e_involved = findall(x->x==s_id, r.reactant_species)
                    if e_involved!=Int64[]
                        push!(s.reaction_list, r)
                    end
                else
                    # For ions and neutral species
                    s_involved = findall(x->x==s_id, r.reactant_species)
                    if s_involved!=Int64[]
                        e_involved = findall(x->x==id_electrons, r.reactant_species)
                        if e_involved==Int64[]
                            # Add only if electron is not involved
                            push!(s.reaction_list, r)
                        end
                    end
                end
            end

            # Flux terms can be considered as a reaction X+ -> X,
            #therefore for each ion (X+) touching the wall a neutral
            #species (X) is created
            if s.charge > 0
                s_neutral = species_list[s.species_id]
                s_neutral.has_wall_loss = true
            end

            # Get pressure for total pressure check
            if system.total_pressure > 0 && s.charge == 0
                pressure_check += s.pressure
            end
        end
        if system.total_pressure > 0
            diff = abs(pressure_check - system.total_pressure)
            if diff > 1.e-10
                print("***ERROR*** Sum of species pressures is not equal to system-defined total pressure\n")
                errcode = c_io_error
            end
        end

    end
    return errcode
end

end