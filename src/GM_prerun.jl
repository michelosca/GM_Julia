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

module GM_prerun

using SharedData: c_io_error, Species, Reaction, System, SpeciesID
using InputData_PreRun: SetupInputData!

function prerun_GM(input_file)

    # Reads data from the input.deck file
    species_list = Species[]
    reaction_list = Reaction[]
    system = System()
    speciesID = SpeciesID()
    errcode = SetupInputData!(input_file, species_list, reaction_list, system,
        speciesID)
    if (errcode == c_io_error)
        return errcode
    end

end

end