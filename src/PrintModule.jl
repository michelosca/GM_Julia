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

module PrintModule

using SharedData: System, Species, Reaction, SpeciesID
using SharedData: e, K_to_eV
using SharedData: h_classical, h_Gudmundsson, h_Monahan 
using SharedData: s_ohmic_power, s_flux_balance, s_flux_interpolation
using InputBlock_Reactions: r_wall_loss, r_elastic
using Printf

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# - PrintSpeciesList
# - PrintReactionList
# - PrintSystemList


function PrintSpeciesList(species_list::Vector{Species}, system::System)

    open(system.log_file,"a") do file
        @printf(file,"Loaded species\n")
        @printf(file, "%15s %15s %15s %15s %10s %10s %10s %10s %15s %15s %15s\n","Name",
            "Species", "Mass [kg]", "Charge [C]","n-eq.",
            "T-eq.", "WL", "P-input","Dens0 [m^-3]","Temp0 [eV]","Press [Pa]")
        for s in species_list

            # s.neutral_id -> species name
            sn_name = species_list[s.species_id].name

            # Density equations
            if (s.has_dens_eq)
                has_dens_eq = "Yes"
            else
                has_dens_eq = "No"
            end

            # Temperature equations
            if (s.has_temp_eq)
                has_temp_eq = "Yes"
            else
                has_temp_eq = "No"
            end

            # Wall losses
            if (s.has_wall_loss)
                has_wall_loss = "Yes"
            else
                has_wall_loss = "No"
            end

            # Power input
            if (s.has_heating_mechanism)
                has_heating = "Yes"
            else
                has_heating = "No"
            end

            @printf(file, "%15s %15s %15.5e %15.5e %10s %10s %10s %10s %15g %15g %15g\n",
                s.name, sn_name, s.mass, s.charge, 
                has_dens_eq, has_temp_eq, has_wall_loss, has_heating,
                s.dens, s.temp*K_to_eV, s.pressure)
        end
        @printf(file, "\n")
    end
end



function PrintReactionList(reaction_list::Vector{Reaction},
    species_list::Vector{Species}, system::System)

    open(system.log_file,"a") do file
        @printf(file, "Loaded reactions\n")
        @printf(file, "%17s %30s %15s %15s %30s %15s\n",
            "ID - description", "Reaction", "E-threshold [eV]",
            "Neutral","Species balance", "Reactants")

        for r in reaction_list

            r_name = string(r.id)
            # r.case -> reaction name
            if (r.case == r_elastic)
                r_name = string(r_name, "-Elastic")
            elseif (r.case == r_wall_loss)
                r_name = string(r_name, "-Wall react.")
            end

            # reaction description
            reactants = String[]
            products = String[]
            for i in 1:length(r.involved_species) 
                s = r.involved_species[i]
                b = r.species_balance[i]

                reactant_species = false 
                s_index = findall( x -> x == s, r.reactant_species )
                if s_index != Int64[] && b > 0 
                    reactant_species = true
                end

                s_name = species_list[s].name
                b_prod = 0
                b_reac = 0
                if b > 0
                    b_prod = b
                end

                if b < 0
                    b_reac = -b
                end

                if reactant_species || b==0
                    b_prod += 1
                    b_reac += 1
                end

                if b_prod == 1
                    push!(products, s_name)
                elseif b_prod > 1
                    push!(products, string(b_prod, s_name))
                end

                if b_reac == 1
                    push!(reactants, s_name)
                elseif b_reac > 1
                    push!(reactants, string(b_reac, s_name))
                end
            end

            # Threshold energy
            E_eV = r.E_threshold / e

            # Neutral species
            r_neutral = ""
            for id in r.neutral_species_id
                r_neutral = string(r_neutral, species_list[id].name,", ")
            end
            r_neutral = chop(r_neutral, tail=2)

            # Involved species and its balance
            r_involved = ""
            i = 1
            for id in r.involved_species
                r_involved = string(r_involved, species_list[id].name,":",r.species_balance[i],", ")
                i += 1
            end
            r_involved = chop(r_involved, tail=2)

            # Reacting species
            r_reactants= ""
            for id in r.reactant_species
                r_reactants = string(r_reactants, species_list[id].name,", ")
            end
            r_reactants = chop(r_reactants, tail=2)

            @printf(file, "%17s %30s %15.2f %15s %30s %15s\n", r_name, r.name,
                E_eV, r_neutral, r_involved, r_reactants)
        end
        @printf(file, "\n")
    end
end


function PrintSystemList(s::System)

    open(s.log_file,"a") do file
        @printf(file, "System parameters\n")
        @printf(file, " - Area:               %15g m2\n", s.A)
        @printf(file, " - Volume:             %15g m3\n", s.V)
        @printf(file, " - Length:             %15g m\n", s.l)
        @printf(file, " - Radius:             %15g m\n", s.radius)
        if (s.h_id == h_classical)
            h_str = "classical"
        elseif (s.h_id == h_Gudmundsson)
            h_str = "Gudmundsson"
        elseif (s.h_id == h_Monahan)
            h_str = "Monahan"
        else
            h_str = "Not defined"
        end
        @printf(file, " - h factor: %s\n",h_str)
        if (s.Vsheath_solving_method== s_ohmic_power)
            vsheath_str = "Ohmic power"
        elseif (s.Vsheath_solving_method== s_flux_balance)
            vsheath_str = "Flux balance (electrons and only negative charged species)"
        elseif (s.Vsheath_solving_method== s_flux_interpolation)
            vsheath_str = "Flux balance (electrons + negative ions)"
        else
            vsheath_str = "Not defined"
        end
        @printf(file, " - Potential solving method: %s\n",vsheath_str)
        @printf(file, " - Driving frequency:      %15g MHz\n", s.drivf/1.e6)
        @printf(file, "                           %15g rad/s\n", s.drivOmega)
        @printf(file, " - Driving power:          %15.2f W\n", s.drivP)
        @printf(file, " - Power shape:            %s \n", s.P_shape)
        @printf(file, " - Power duty ratio:       %15g \n", s.P_duty_ratio)
        @printf(file, " - Total neutral pressure: %15g Pa\n", s.total_pressure)
        @printf(file, " - Plasma potential guess: %15g V\n", s.plasma_potential)
        @printf(file, " - Simulation time:        %15g s\n", s.t_end)
        @printf(file, " - Simulation folder:      %s\n", s.folder)
        @printf(file, " - Log file:               %s\n", s.log_file)
        @printf(file, "\n")
    end
end

function PrintSimulationState(temp::Vector{Float64}, dens::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)
    PrintSystemList(system)
    PrintSpeciesList(species_list, system)
    open(system.log_file,"a") do file
        @printf(file,"Plasma potential: %15g\n", system.plasma_potential) 
        @printf(file,"Total pressure:   %15g\n", system.total_pressure) 
        @printf(file,"Electronegativity:%15g\n", system.alpha) 
        @printf(file,"Species values\n")
        for s in species_list
            @printf(file,"  - Species %s\n    - array-dens: %15g\n    - array-temp: %15g\n", s.name, dens[s.id], temp[s.id])
            @printf(file,"    - pressure: %15g\n",s.pressure) 
            @printf(file,"    - flux: %15g\n",s.flux) 
            @printf(file,"    - mfp: %15g\n", s.mfp) 
        end
    end
end

function PrintErrorMessage(system::System, message::String)

    open(system.log_file,"a") do file
        @printf(file, "***ERROR*** %s\n", message)
    end
    @printf("***ERROR*** %s\n", message)
end

function PrintWarningMessage(system::System, message::String)

    open(system.log_file,"a") do file
        @printf(file, "***WARNING*** %s\n", message)
    end
    @printf("***WARNING*** %s\n", message)
end

end