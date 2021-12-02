module PlotModule

using Plots
#using LaTeXStrings
using SharedData: Species, OutputBlock, Reaction
using SharedData: K_to_eV

function PlotDensities_Charged(output::OutputBlock,
    species_list::Vector{Species})

    p = plot(xlabel = "t [ms]", ylabel = "n [m^-3]",
        xscale = :log10, xlims = (5.e-2, 1.e4), yscale = :log10,
        ylims = (1.e14, 1.e17))

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

    time = output.x/1.e-3
    for s in species_list
        if !s.has_dens_eq
            continue
        end
        if s.charge == 0
            dens = output.n[s.id]
            plot!(p, time, dens, label = s.name, w = 2)
        end
    end
    return p
end


function PlotTemperatures_Charged(output::OutputBlock,
    species_list::Vector{Species})

    p = plot(xlabel = "t [ms]", ylabel = "T [eV]",
        xscale = :log10, xlims = (5.e-2, 1.e4),
        ylims = (0, 5))

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

    p = plot(xlabel = "t [ms]", ylabel = "m^3/s",
        yscale = :log10, ylims = (1.e-20, 1.e-10), legend=:outertopright)#,
        #ylims = (0, 5))
    time = output.x#/1.e-3
    if reaction_id == 0
        for r in reaction_list 
            K = output.K[r.id]
            if K == Float64[]
                continue
            else
                plot!(p, time, K, label=r.name, w = 2, yscale=:log10)
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