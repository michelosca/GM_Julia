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

module PlotModule

using Plots
using Plots.PlotMeasures
#using LaTeXStrings
using SharedData: Species, OutputBlock, Reaction
using SharedData: K_to_eV
using DataFrames

function Get2DMap(df::DataFrame, name::String)

    # First set M_col and M_row for matrix M
    dat = df[1,2]
    row_flag = true
    M_col = 0
    M_row = 0
    for row in 1:nrow(df)
        if dat != df[row,2]
            if row_flag
                row_flag = false
                M_row = row - 1
            end
            M_col += 1
            dat = df[row,2]
        end
    end
    M_col += 1

    # Fill matrix M and col/row vectors
    M = zeros(M_row, M_col)
    row_vec = zeros(M_row)
    col_vec = zeros(M_col)
    dat = df[1,2]
    col = 1
    col_vec[col] = dat
    row_flag = true
    for row in 1:nrow(df)
        # Fill coll vector
        if dat != df[row,2]
            row_flag = false
            col += 1
            dat = df[row,2]
            col_vec[col] = dat
        end
        
        # Fill row vector
        if row_flag
            row_vec[row] = df[row, 1]
        end

        # Fill matrix
        current_row = row - M_row*(col-1)
        current_col = col
        M[current_row, current_col] = df[row,name]
    end

    return row_vec, col_vec, M
end

function PlotHeatmap(row, col, matrix, title_str)

    h = heatmap(row,
    col, matrix,
    size=(1400, 800),
    ylabel="O2 fraction",
    xlabel = "Power [W]",
    title = title_str,
    #colorbar_title="m^3/s",
    xtickfont=20,
    ytickfont=20,
    titlefont=30,
    yguidefontsize=20,
    xguidefontsize=20,
    left_margin=10mm,
    right_margin=40mm,
    bottom_margin=10mm,
    top_margin=10mm,
    colorbar_titlefontsize=20)#colorbar_title_locations=20mm)

    return h
end


function PlotDensities_Charged(output::OutputBlock,
    species_list::Vector{Species})

    p = plot()#xlabel = "t [ms]",
        #ylabel = "n [m^-3]",
        #xscale = :log10,
        #xlims = (5.e-2, 1.e4),
        #yscale = :log10,
        #ylims = (1.e14, 1.e17))

    time = output.x#/1.e-3
    for s in species_list
        if !s.has_dens_eq
            continue
        end
        if abs(s.charge) > 0
            dens = output.n[s.id]
            plot!(p, time, dens, label=s.name, w = 2)
        end
    end
    return p
end


function PlotDensities_Neutral(output::OutputBlock,
    species_list::Vector{Species})

    p = plot(xlabel = "t [ms]", ylabel = "n [m^-3]")#,
    #yscale = :log10, xlims = (1.e15, 1.e22))

    time = output.n_data_frame[!,1]
    for s in species_list
        if !s.has_dens_eq
            continue
        end
        if s.charge == 0
            dens = output.n_data_frame[!,s.name]
            plot!(p, time, dens, label = s.name, w = 2)
        end
    end
    return p
end


function PlotTemperatures_Charged(output::OutputBlock,
    species_list::Vector{Species})

    p = plot(xlabel = "t [ms]", ylabel = "T [eV]")#,
        #xscale = :log10, xlims = (5.e-2, 1.e4),
        #ylims = (0, 5))

    time = output.x#/1.e-3
    for s in species_list
        if !s.has_temp_eq
            continue
        end
        if abs(s.charge) > 0
            temp = output.T[s.id] * K_to_eV
            plot!(p, time, temp, label=s.name, w = 2)
        end
    end
    return p
end


function PlotTemperatures_Neutral(output::OutputBlock,
    species_list::Vector{Species})

    p = plot(xlabel = "t [ms]", ylabel = "T [eV]")

    time = output.x/1.e-3
    for s in species_list
        if !s.has_temp_eq
            continue
        end
        if s.charge == 0
            temp = output.n[s.id] * K_to_eV
            plot!(p, time, temp, label = s.name, w = 2)
        end
    end
    return p
end


function PlotRateCoefficients(output::OutputBlock,
    reaction_list::Vector{Reaction}, reaction_id::Int64 = 0)

    p = plot(xlabel = "Reaction ID", ylabel = "",
        kind = "bar",
        yscale = :log10)
        #ylims = (1.e-20, 1.e-10),
        #legend=:outertopright,
        #ylims = (0, 5))
    time = output.x#/1.e-3
    if reaction_id == 0
        for r in reaction_list 
            K = output.K[r.id]
            if K == Float64[]
                continue
            elseif K[end] == 0
                continue
            else
                dens_list = output.n[r.reactant_species]
                dens = 1
                for n in dens_list
                    dens *= n[end]
                end
                #print("ID ", r.id, " dens product ", dens,"  K ", K[end]," Kn ",K[end]*dens ,"\n")
                plot!(p, [r.id], [K[end]*dens], w = 2, markershape=:circle, legend = false)
            end
        end
    else
        r = reaction_list[reaction_id]
        K = output.K[r.id]
        plot!(p, time, K, label=r.name, w = 2, yscale=:log10)
    end

    return p
end

end