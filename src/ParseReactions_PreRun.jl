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

module ParseReactions_PreRun

using SharedData: c_io_error
using SharedData: r_extended, r_diffusion, r_emission_rate
using SharedData: Reaction

using EvaluateExpressions: ReplaceConstantValues!, ReplaceSystemSymbols!
using EvaluateExpressions: ReplaceSpeciesSymbols!, ReplaceTempSymbols!
using EvaluateExpressions: ReplaceDensSymbols!



function WriteRateCoefficientFunctions!(current_reaction::Reaction,
    reaction_str::SubString{String}, f_ReactionSet::IOStream)

    errcode = 0
    
    try
        expr = Meta.parse(reaction_str)
        if !(typeof(expr)==Float64)
            ReplaceConstantValues!(expr)
            ReplaceSystemSymbols!(expr)
            ReplaceSpeciesSymbols!(expr)
            ReplaceTempSymbols!(expr)
            ReplaceDensSymbols!(expr)
        end

        # Now that the expression is ready to be evaluated, write it down in a new file
        if current_reaction.case == r_extended || 
            current_reaction.case == r_diffusion ||
            current_reaction.case == r_emission_rate
            write(f_ReactionSet,
                string("push!(K_funct_list, (dens::Vector{Float64}, ",
                "temp::Vector{Float64}, species_list::Vector{Species}, ",
                "system::System,  sID::SpeciesID) -> ", expr,")\n"))
        else
            write(f_ReactionSet,
                string("push!(K_funct_list, (temp::Vector{Float64},",
                " sID::SpeciesID) -> ",expr,")\n"))
        end

    catch
        errcode = c_io_error
        print("***ERROR*** While writing rate coefficient to module\n")
    end

    return errcode
end

end