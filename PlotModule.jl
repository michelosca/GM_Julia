module PlotModule

using Plots
#using LaTeXStrings

function PlotData(x,y)
    p = plot(x, y,
    xlabel = "p·L [Pa·m]",
    #ylabel = "Electron Density [m^{-3}]",
    ylabel = "Electron Temp [eV]",
    title = "Ar discharge 1kW @ 200V",
    xlims = (5.e-2,1.e4),
    #ylims = (1.e14, 1.e17),
    ylims = (0, 5),
    w = 4,
    grid = false,
    xticks = [1.e-1,1.e0,1.e1,1.e2,1.e3,1.e4],
    #yticks = [1.e14,1.e15,1.e16,1.e17],
    xscale = :log10,
    #yscale = :log10,
    #xtickfont = font(20, "Courier"),
    #ytickfont = font(20, "Courier"),
    legend = false,
    framestyle = :box
    )
    return p
end

end