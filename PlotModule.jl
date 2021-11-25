module PlotModule

using Plots
#using LaTeXStrings
using SharedData: Species

function PlotDens(x,y)
    p = plot(x, y,
    xlabel = "p·L [Pa·m]",
    ylabel = "Electron Density [m^{-3}]",
    title = "Ar discharge 1kW @ 200V",
    xlims = (5.e-2,1.e4),
    ylims = (1.e14, 1.e17),
    w = 4,
    linestyle = :dot,
    grid = false,
    xticks = [1.e-1,1.e0,1.e1,1.e2,1.e3,1.e4],
    yticks = [1.e14,1.e15,1.e16,1.e17],
    xscale = :log10,
    minorticks = true,
    yscale = :log10,
    #xtickfont = font(20, "Courier"),
    #ytickfont = font(20, "Courier"),
    legend = false,
    framestyle = :box
    )
    return p
end

function PlotTemp(x,y)
    p = plot(x, y,
    xlabel = "p·L [Pa·m]",
    ylabel = "Electron Temp [eV]",
    title = "Ar discharge 1kW @ 200V",
    xlims = (5.e-2,1.e4),
    ylims = (0, 5),
    w = 4,
    linestyle = :dot,
    grid = false,
    xticks = [1.e-1,1.e0,1.e1,1.e2,1.e3,1.e4],
    yticks = [0,1,2,3,4,5],
    xscale = :log10,
    minorticks = true,
    #xtickfont = font(20, "Courier"),
    #ytickfont = font(20, "Courier"),
    legend = false,
    framestyle = :box
    )
    return p
end

function PlotDensities_Charged(species_list::Vector{Species}, sol)
    # Input: sol, is the output returned by DifferentialEquations.solve
    p = plot(xlabel = "t [ms]", ylabel = "n [m^-3]")
    n_species = length(species_list)
    for s in species_list
        if abs(s.charge) > 0
            dens_id = s.id + n_species
            dens = sol[dens_id, :]
            time = sol.t/1.e-3
            plot!(p, time, dens, label=s.name, w = 2)
        end
    end
    return p
end

function PlotDensities_Neutral(species_list::Vector{Species}, sol)
    # Input: sol, is the output returned by DifferentialEquations.solve
    p = plot(xlabel = "t [ms]", ylabel = "n [m^-3]")
    n_species = length(species_list)
    for s in species_list
        if s.charge == 0
            dens_id = s.id + n_species
            dens = sol[dens_id, :]
            time = sol.t/1.e-3
            plot!(p, time, dens, label = s.name, w = 2)
        end
    end
    return p
end

end