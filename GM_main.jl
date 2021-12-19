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

module GM_main

using SharedData: c_io_error, Species, Reaction, System, SpeciesID, OutputBlock
using InputData: SetupInputData!
using PrintModule: PrintSpeciesList, PrintReactionList, PrintSystemList
using OutputModule: GenerateOutputs!
using CSV

function run_GM(input)

    if typeof(input) == String
        # Reads data from the input.deck file
        species_list = Species[]
        reaction_list = Reaction[]
        system = System()
        output_list = OutputBlock[]
        speciesID = SpeciesID()
        errcode = SetupInputData!(input, species_list,
            reaction_list, system, output_list, speciesID)
        if (errcode == c_io_error)
            return species_list, reaction_list, system, speciesID, output_list 
        end

    elseif typeof(input) == Tuple{Vector{Species}, Vector{Reaction}, System,
        Vector{OutputBlock}, SpeciesID} 
        # Restart from previous simulation
        species_list = input[1]
        reaction_list = input[2]
        system = input[3]
        output_list = input[4]
        speciesID = input[5]
    else
        print("***ERROR*** Input is not recognized\n")
        if (errcode == c_io_error)
            return nothing, nothing, nothing, nothing, nothing
        end
    end

    # Print system, species and reaction lists to terminal
    PrintSystemList(system)
    PrintSpeciesList(species_list, speciesID)
    PrintReactionList(reaction_list, species_list, speciesID)

    errcode = @time GenerateOutputs!(species_list, reaction_list, system,
        output_list, speciesID)
    if (errcode == c_io_error)
        return species_list, reaction_list, system, speciesID, output_list 
    end

    for output in output_list
        label = ""
        for i in 1:output.n_parameters
            label = string(label,"_vs_",output.name[i])
        end
        T_filename = string("T",label,".csv")
        n_filename = string("n",label,".csv")
        K_filename = string("K",label,".csv")
        CSV.write(T_filename, output.T_data_frame )
        CSV.write(n_filename, output.n_data_frame )
        CSV.write(K_filename, output.K_data_frame )
    end

    return species_list, reaction_list, system, speciesID, output_list
end

end