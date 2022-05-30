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
using SharedData: Species, OutputBlock, Reaction, SpeciesID, System
using SharedData: K_to_eV, kb
using DataFrames

function CRM_test(out_batch::Vector{Vector{OutputBlock}})

    p_total = 5.0

    # %O2: 2.5/102.5, 5/105, 7.5/107.5,10/110, 20/120
    O2_fract = [2.5/102.5, 5/105, 7.5/107.5,10/110, 20/120, 1]
    T_g_fieb2020_K = [622,683,717,739,787,700]
    Te_fieb2020_eV = [1.9, 2.2, 2.5, 2.6, 2.7, 2]
    ne_fieb2020_e16 = [52.10, 33.00, 21.23, 16.29, 10.84, 5.18]
    Ar_m_fieb2020_e16 = [2.2+0.23, 2.8+0.28, 2.4+0.16, 1.9+0.08, 0.31+0.01, 0.0]
    Ar_r_fieb2020_e16 = [0.68+0.55, 0.82+0.67, 0.61+0.47, 0.46+0.32, 0.07+0.05, 0.0]
    I_130_fieb2020_e18 = [450, 510, 460, 400, 260,82]
    I_135_fieb2020_e18 = [29,37,40,35,18,1.3]
    I_777_fieb2020_e18 = [520,570,470,340,150,32]
    nO2_fieb2020_e20 = [0.1, 0.17, 0.25, 0.33, 0.64, 5.02]
    nO_fieb2020_e19 = [0.41, 0.93, 1.4, 1.7, 2.3, 3.15]
    nOall_fieb2020_e19 = [0.83, 1.55, 2.0, 2.23, 2.62, 3.17]
    nO_1d_fieb2020_e18 = [4.1, 6.1, 5.9, 5.2, 3.1, 0.23]
    nO_1s_fieb2020_e17 = [1.3, 1.5, 0.79, 0.73, 0.54,0.11]
    nO_5s_fieb2020_e15 = [5.2, 6.7, 7.2, 6.3, 3.2, 0.23]
    nO_3s_fieb2020_e13 = [0.6, 1.4, 1.7, 1.8, 1.5, 0.66]
    nO_5p_fieb2020_e13 = [2.9, 3.5, 3.0, 2.0, 0.65, 0.09]
    nO_3p_fieb2020_e12 = [4.7, 6.9, 9.0, 9.2, 7.2, 1.1]
    dissociation_fieb2020_percent = [29.5, 30.8, 28.3, 25.1, 17.1, 3.1]

    Te_list = Vector{Float64}[]
    ne_list = Vector{Float64}[]
    Ar_m_list = Vector{Float64}[]
    Ar_r_list = Vector{Float64}[]
    I_130_list = Vector{Float64}[]
    I_135_list = Vector{Float64}[]
    I_777_list = Vector{Float64}[]
    nO2_list = Vector{Float64}[]
    nO_list = Vector{Float64}[]
    nO_1d_list = Vector{Float64}[]
    nO_1s_list = Vector{Float64}[]
    nO_5s_list = Vector{Float64}[]
    nO_3s_list = Vector{Float64}[]
    nO_5p_list = Vector{Float64}[]
    nO_3p_list = Vector{Float64}[]
    dissociation_list = Vector{Float64}[] 

    for out_list in out_batch
        Te_GM_eV = Float64[]
        ne_GM_e16 = Float64[]
        Ar_m_GM_e16 = Float64[]
        Ar_r_GM_e16 = Float64[]
        I_130_GM_e18 = Float64[]
        I_135_GM_e18 = Float64[]
        I_777_GM_e18 = Float64[]
        nO2_GM_e20 = Float64[]
        nO_GM_e19 = Float64[]
        nO_1d_GM_e18 = Float64[]
        nO_1s_GM_e17 = Float64[]
        nO_5s_GM_e15 = Float64[]
        nO_3s_GM_e13 = Float64[]
        nO_5p_GM_e13 = Float64[]
        nO_3p_GM_e12 = Float64[]
        dissociation_GM_percent = Float64[] 

        for i in range(1, length(out_list), step=1)
            out = out_list[i]
            n_data = out.n_data_frame
            T_data = out.T_data_frame
            K_data = out.K_data_frame
            nO_all = 0.0

            #print("    T_e: ", T_data.e[end]*K_to_eV,"\n")
            push!(Te_GM_eV, T_data.e[end]*K_to_eV)

            #print("    n_e: ", n_data.e[end],"\n")
            push!(ne_GM_e16, n_data.e[end]*1.e-16)

            #print("   n_02: ", n_O2,"\n")
            n_O2 = n_data.O2[end]
            push!(nO2_GM_e20, n_O2*1.e-20) 

            nO = n_data.O[end]
            nO_all += nO
            #print("    n_0: ", nO,"\n")
            push!(nO_GM_e19, nO*1.e-19)

            #print(" n_0_1d: ", n_data.O_1d[end],"\n")
            nO_all += n_data.O_1d[end]
            push!(nO_1d_GM_e18, n_data.O_1d[end]*1.e-18)

            #print(" n_0_1s: ", n_data.O_1s[end],"\n")
            nO_all += n_data.O_1s[end]
            push!(nO_1s_GM_e17, n_data.O_1s[end]*1.e-17)

            #print(" n_0_5s: ", n_data.O_5s[end],"\n")
            nO_all += n_data.O_5s[end]
            push!(nO_5s_GM_e15, n_data.O_5s[end]*1.e-15)

            #print(" n_0_3s: ", n_data.O_3s[end],"\n")
            nO_all += n_data.O_3s[end]
            push!(nO_3s_GM_e13, n_data.O_3s[end]*1.e-13)

            #print(" n_0_5p: ", n_data.O_5p[end],"\n")
            nO_all += n_data.O_5p[end]
            push!(nO_5p_GM_e13, n_data.O_5p[end]*1.e-13)

            #print(" n_0_3p: ", n_data.O_3p[end],"\n")
            nO_all += n_data.O_3p[end]
            push!(nO_3p_GM_e12, n_data.O_3p[end]*1.e-12)
            
            #print("  I_130: ", K_data."r300: O_3s -> O"[end],"\n")
            I_130nm = K_data."r297: O_3s -> O"[end]
            I_130nm += K_data."r332: e + O2 -> e + 2O"[end]
            #I_130nm += K_data."r333: e + O2_v -> e + 2O"[end]
            I_130nm += K_data."r333: e + O2_a1Ag -> e + 2O"[end]
            I_130nm += K_data."r334: e + O2_b1Su -> e + 2O"[end]
            I_130nm += K_data."r335: e + O2_a1Ag_v -> e + 2O"[end]
            I_130nm += K_data."r336: e + O2_b1Su_v -> e + 2O"[end]
            push!(I_130_GM_e18, I_130nm*1.e-18)

            #print("  I_135: ", K_data."r299: O_5s -> O"[end],"\n")
            I_135nm = K_data."r296: O_5s -> O"[end]
            I_135nm += K_data."r337: e + O2 -> e + 2O"[end]
            #I_135nm += K_data."r338: e + O2_v -> e + 2O"[end]
            I_135nm += K_data."r338: e + O2_a1Ag -> e + 2O"[end]
            I_135nm += K_data."r339: e + O2_b1Su -> e + 2O"[end]
            I_135nm += K_data."r340: e + O2_a1Ag_v -> e + 2O"[end]
            I_135nm += K_data."r341: e + O2_b1Su_v -> e + 2O"[end]
            push!(I_135_GM_e18, I_135nm*1.e-18)
            
            #print("  I_777: ", K_data."r301: O_5p -> O_5s"[end],"\n")
            I_777nm = K_data."r298: O_5p -> O_5s"[end]
            I_777nm += K_data."r331: e + O -> e + O_5s"[end]
            I_777nm += K_data."r342: e + O2 -> e + O + O_5s"[end]
            I_777nm += K_data."r343: e + O2_a1Ag -> e + O + O_5s"[end]
            I_777nm += K_data."r344: e + O2_b1Su -> e + O + O_5s"[end]
            I_777nm += K_data."r345: e + O2_a1Ag_v -> e + O + O_5s"[end]
            I_777nm += K_data."r346: e + O2_b1Su_v -> e + O + O_5s"[end]
            push!(I_777_GM_e18, I_777nm*1.e-18)

            #print("Ar_m 5s: ", n_data.Ar_m[end],"\n")
            push!(Ar_m_GM_e16, n_data.Ar_m[end]*1.e-16)

            #print("Ar_r 4s: ", n_data.Ar_r[end],"\n")
            push!(Ar_r_GM_e16, n_data.Ar_r[end]*1.e-16)
            
            #print("Dissoci: ", 100* 0.5*nO / (n_O2 + 0.5*nO),"\n")
            push!(dissociation_GM_percent, 100* 0.5*nO_all / (n_O2 + 0.5*nO_all))
        end

        push!(Te_list, Te_GM_eV)
        push!(ne_list, ne_GM_e16)
        push!(nO2_list, nO2_GM_e20)
        push!(nO_list, nO_GM_e19)
        push!(nO_1d_list, nO_1d_GM_e18)
        push!(nO_1s_list, nO_1s_GM_e17)
        push!(nO_3s_list, nO_3s_GM_e13)
        push!(nO_5s_list, nO_5s_GM_e15)
        push!(nO_3p_list, nO_3p_GM_e12)
        push!(nO_5p_list, nO_5p_GM_e13)
        push!(Ar_m_list, Ar_m_GM_e16)
        push!(Ar_r_list, Ar_r_GM_e16)
        push!(I_130_list, I_130_GM_e18)
        push!(I_135_list, I_135_GM_e18)
        push!(I_777_list, I_777_GM_e18)
        push!(dissociation_list, dissociation_GM_percent)
    end

    source = "/home/moe505/Documents/GM_Julia/ArO2/solve_O2/"
    #str_list = ["gamma_O2=0.0", "gamma_O2=0.007", "gamma_02=0.1", "gamma_02=1.0"]
    #str_list = ["gamma_02=0.007 / O_x -> O", "gamma_02=0.007 / 2O_x -> O2"]
    #str_list = ["gamma_Ar=1.0", "gamma_Ar=0.0"]
    #str_list = ["gamma_O-meta=0.0", "gamma_O-meta=0.0 / 2O_x -> O2", "gamma_O-meta=1.0", "gamma_O-meta=1.0 / 2O_x -> O2"]
    str_list = ["GM Julia","GM Julia, solve O2"]
    color_list = [:black, :blue, :green, :red]

    # ELECTRON TEMPERATURE
    p_Te = scatter(O2_fract, Te_fieb2020_eV, label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_Te, O2_fract, Te_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_Te, title="electron temperature", xlabel="O2 fraction", ylabel="eV",
        grid=true, xscale=:log10)#, ylim=(1.5,3))
    png(source*"Te_vs_Fieb2020")

    # O2 DENSITY
    p_nO2 = scatter(O2_fract, nO2_fieb2020_e20, label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_nO2, O2_fract, nO2_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_nO2, title="O2 density", xlabel="O2 fraction", ylabel="1e20/m^3",
        grid=true, xscale=:log10, legend=:topleft, yscale=:log10, ylim=(0.01,10)) 
    png(source*"nO2_vs_Fieb2020")
    
    # O2-DISSOCIATION RATE
    p_alpha = scatter(O2_fract, dissociation_fieb2020_percent, label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_alpha, O2_fract, dissociation_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_alpha, title="O/O2 dissociation rate [%]", xlabel="O2 fraction", ylabel="%",
        grid=true, xscale=:log10, legend=:topright, yscale=:identity, ylim=(0,40)) 
    png(source*"dissociation_vs_Fieb2020")
    
    # O DENSITY
    p_nO = scatter(O2_fract, nO_fieb2020_e19, label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_nO, O2_fract, nO_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_nO, title="O density", xlabel="O2 fraction", ylabel="1e19/m^3",
        grid=true, xscale=:log10, legend=:bottomright, yscale=:log10, ylim=(0.1,10))
    png(source*"nO_vs_Fieb2020")
    
    # O_1d DENSITY
    #p_nO1d = scatter(O2_fract, nO_1d_fieb2020_e18, label="Fiebrandt 2020",
    #    markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    #for i in range(1, length(out_batch), step=1)
    #    scatter!(p_nO1d, O2_fract, nO_1d_list[i], label=str_list[i],
    #        markercolor=i, markerstrokecolor=i, markersize=6)
    #end
    #scatter!(p_nO1d, title="O 1d density", xlabel="O2 fraction", ylabel="1e18/m^3",
    #    grid=true, xscale=:log10, legend=:bottomright, ylim=(0.1,10), yscale=:log10)
    #png(source*"nO1d_vs_Fieb2020")
    
    # O_1s DENSITY
    #p_nO1s = scatter(O2_fract, nO_1s_fieb2020_e17, label="Fiebrandt 2020",
    #    markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    #for i in range(1, length(out_batch), step=1)
    #    scatter!(p_nO1s, O2_fract, nO_1s_list[i], label=str_list[i],
    #        markercolor=i, markerstrokecolor=i, markersize=6)
    #end
    #scatter!(p_nO1s, title="O 1s density", xlabel="O2 fraction", ylabel="1e17/m^3",
    #    grid=true, xscale=:log10, legend=:bottomright, ylim=(0.01,10), yscale=:log10)
    #png(source*"nO1s_vs_Fieb2020")
    
    # O_5s DENSITY
    #p_nO5s = scatter(O2_fract, nO_5s_fieb2020_e15, label="Fiebrandt 2020",
    #    markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    #for i in range(1, length(out_batch), step=1)
    #    scatter!(p_nO5s, O2_fract, nO_5s_list[i], label=str_list[i],
    #        markercolor=i, markerstrokecolor=i, markersize=6)
    #end
    #scatter!(p_nO5s, title="O 5s density", xlabel="O2 fraction", ylabel="1e15/m^3",
    #    grid=true, xscale=:log10, legend=:bottomright, ylim=(0.1,10), yscale=:log10)
    #png(source*"nO5s_vs_Fieb2020")
    
    #p_nO3s = scatter(O2_fract, nO_3s_fieb2020_e13, label="Fiebrandt 2020",
    #    markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    #for i in range(1, length(out_batch), step=1)
    #    scatter!(p_nO3s, O2_fract, nO_3s_list[i], label=str_list[i],
    #        markercolor=i, markerstrokecolor=i, markersize=6)
    #end
    #scatter!(p_nO3s, title="O 3s density", xlabel="O2 fraction", ylabel="1e13/m^3",
    #    grid=true, xscale=:log10, legend=:bottomright, ylim=(0.01,10), yscale=:log10)
    #png(source*"nO3s_vs_Fieb2020")
    
    #p_nO5p = scatter(O2_fract, nO_5p_fieb2020_e13, label="Fiebrandt 2020",
    #    markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    #for i in range(1, length(out_batch), step=1)
    #    scatter!(p_nO5p, O2_fract, nO_5p_list[i], label=str_list[i],
    #        markercolor=i, markerstrokecolor=i, markersize=6)
    #end
    #scatter!(p_nO5p, title="O 5p density", xlabel="O2 fraction", ylabel="1e13/m^3",
    #    grid=true, xscale=:log10, legend=:bottomright, ylim=(0.01,10), yscale=:log10)
    #png(source*"nO5p_vs_Fieb2020")
    
    #p_nO3p = scatter(O2_fract, nO_3p_fieb2020_e12, label="Fiebrandt 2020",
    #    markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    #for i in range(1, length(out_batch), step=1)
    #    scatter!(p_nO3p, O2_fract, nO_3p_list[i], label=str_list[i],
    #        markercolor=i, markerstrokecolor=i, markersize=6)
    #end
    #scatter!(p_nO3p, title="O 3p density", xlabel="O2 fraction", ylabel="1e12/m^3",
    #    grid=true, xscale=:log10, legend=:bottomright, ylim=(0.1,10), yscale=:log10)
    #png(source*"nO3p_vs_Fieb2020")
    
    # Ar_m DENSITY
    p_Arm = scatter(O2_fract[1:5], Ar_m_fieb2020_e16[1:5], label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_Arm, O2_fract[1:5], Ar_m_list[i][1:5], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_Arm, title="Ar-m density", xlabel="O2 fraction", ylabel="1e16/m^3",
        grid=true, xscale=:identity, legend=:topright, yscale=:log10, ylim=(0.1,10))
    png(source*"Ar_m_vs_Fieb2020")
    
    # Ar_r DENSITY
    p_Arr = scatter(O2_fract[1:5], Ar_r_fieb2020_e16[1:5], label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_Arr, O2_fract[1:5], Ar_r_list[i][1:5], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_Arr, title="Ar-r density", xlabel="O2 fraction", ylabel="1e16/m^3",
        grid=true, xscale=:identity, legend=:topright, yscale=:log10, ylim=(0.1,10)) 
    png(source*"Ar_r_vs_Fieb2020")
    
    # ELECTRON DENSITY
    p_ne = scatter(O2_fract, ne_fieb2020_e16, label="Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_ne, O2_fract, ne_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_ne, title="electron density", xlabel="O2 fraction", ylabel="1e16/m^3",
        grid=true, xscale=:log10, legend=:topright, yscale=:log10, ylim=(1,100))
    png(source*"ne_vs_Fieb2020")
    
    # EMISSION RATE AT 130nm
    p_I130 = scatter(O2_fract, I_130_fieb2020_e18, label="130nm Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_I130, O2_fract, I_130_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_I130, title="emission 130", xlabel="O2 fraction", ylabel="1e18/s/m^3",
        grid=true, xscale=:log10, ylim=(10,1000), yscale=:log10)
    png(source*"I130nm_vs_Fieb2020")

    # EMISSION RATE AT 135nm
    p_I135 = scatter(O2_fract, I_135_fieb2020_e18, label="135nm Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_I135, O2_fract, I_135_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_I135, title="emission 135", xlabel="O2 fraction", ylabel="1e18/s/m^3",
        grid=true, xscale=:log10, ylim=(1,100), yscale=:log10)
    png(source*"I135nm_vs_Fieb2020")

    # EMISSION RATE AT 777nm
    p_I777 = scatter(O2_fract, I_777_fieb2020_e18, label="777nm Fiebrandt 2020",
        markercolor=:red, markerstrokecolor=:red, markersize=6, marker=:cross)
    for i in range(1, length(out_batch), step=1)
        scatter!(p_I777, O2_fract, I_777_list[i], label=str_list[i],
            markercolor=i, markerstrokecolor=i, markersize=6)
    end
    scatter!(p_I777, title="emission 777", xlabel="O2 fraction", ylabel="1e18/s/m^3",
        grid=true, xscale=:log10, ylim=(10,1000), yscale=:log10)
    png(source*"I777nm_vs_Fieb2020")
end


function FindProductSpecies(reaction_list::Vector{Reaction}, species::Species)
    sid = species.id
    output_reaction_list = Reaction[]
    for r in reaction_list
        # Get species index in involved_species list
        s_index = findall( x -> x == sid, r.involved_species )
        if s_index == Int64[]
            continue
        else
            s_index = s_index[1]
        end

        sign = r.species_balance[s_index]
        if sign > 0
            print(r.name, "\n")
            push!(output_reaction_list, r) 
        end
    end
    return output_reaction_list
end

function FindReactingSpecies(reaction_list::Vector{Reaction}, species::Species)
    sid = species.id
    output_reaction_list = Reaction[]
    for r in reaction_list
        # Get species index in involved_species list
        s_index = findall( x -> x == sid, r.involved_species )
        if s_index == Int64[]
            continue
        else
            s_index = s_index[1]
        end

        sign = r.species_balance[s_index]
        if sign < 0
            print(r.name, "\n")
            push!(output_reaction_list, r) 
        end
    end
    return output_reaction_list
end


function PlotSpeciesReactionFraction(K_df::DataFrame, species::Species, reaction_list::Vector{Reaction})
    source = "/home/moe505/Documents/GM_Julia/ArO2/Param_Investigation:P_total-Power/"

    col1 = 1
    col2 = 3
    # Dummy call to get matrix size
    row, col, mat = Get2DMap(K_df, names(K_df)[1],col1,col2)
    mat_tot = zeros(length(row), length(col))

    # Identify reactions with given species in RHS
    # Add all these reactions to the mat_tot matrix
    i = 0 
    species_reaction_list = FindProductSpecies(reaction_list, species)
    for r in species_reaction_list
        name= names(K_df)[r.id + col2]
        row, col, mat = Get2DMap(K_df, name, col1, col2)
        mat_tot .+= mat
    end

    # Plot reactions
    for r in species_reaction_list
        reaction_str = names(K_df)[r.id + col2]
        row, col, mat = Get2DMap(K_df, reaction_str,col1,col2)
        
        P_tot = row
        Power = col
        PlotHeatmap(P_tot, Power, mat./mat_tot, reaction_str )
        png(source*reaction_str*"_fraction")
    end
end


function Plot2DHeatmap(df::DataFrame, flag::Int64)

    source = "/home/moe505/Documents/GM_Julia/ArO2/Param_Investigation:P_total-Power/"


    col1 = 1 # pressure
    col2 = 3 # power
    if flag == 1  # density data frame
        for s_name in names(df)[col2+1:end] 
            print(s_name,"\n")
            label = "number density (log10 scale): " * s_name
            row, col, matrix = Get2DMap(df, s_name,col1,col2)
            PlotHeatmap(row, col, log10.(matrix), label)
            fig_name = source * "dens_" * s_name
            png(fig_name)
        end
    elseif flag == 2  # rate coefficient data frame
        for s_name in names(df)[col2+1:end] 
            print(s_name,"\n")
            label = "Collision rate [1/s] (log10 scale): " * s_name
            row, col, matrix = Get2DMap(df, s_name,col1,col2)
            PlotHeatmap(row, col, log10.(matrix), label)
            fig_name = source * s_name
            png(fig_name)
        end
    end

end


function Get2DMap(df::DataFrame, name::String, col1::Int64, col2::Int64)

    # First: index the row and column arrays 
    M_row = Float64[df[1,col1]] 
    M_col = Float64[df[1,col2]]
    for df_row in 1:nrow(df)
        df_row_value = df[df_row, col1]
        df_col_value = df[df_row, col2]
        row_is_contained = df_row_value in M_row
        col_is_contained = df_col_value in M_col
        if !row_is_contained
            push!(M_row, df_row_value)
        end
        if !col_is_contained
            push!(M_col, df_col_value)
        end
    end
    M_row = sort(M_row)
    M_col = sort(M_col)

    # Fill matrix M and col/row vectors
    n_rows = length(M_row)
    n_cols = length(M_col)
    M = zeros(n_rows, n_cols)
    for df_row in 1:nrow(df)
        df_row_value = df[df_row, col1]
        df_col_value = df[df_row, col2]
        row = findfirst(x -> x==df_row_value, M_row)
        col = findfirst(x -> x==df_col_value, M_col)
        df_value = df[df_row, name]
        if (M[row, col] == 0.0)
            M[row, col] = df_value 
        else
            M[row, col] += df_value 
            M[row, col] *= 0.5 
        end
    end

    return M_row, M_col, M
end


function PlotHeatmap(row, col, matrix, title_str)

    #h = heatmap(col,
    h = contourf(col,
    row, matrix,
    size=(1400, 800),
    ylabel= "P total [Pa]",     # row data label
    xlabel = "Power [W]",  # row data label
    xscale = :identity,
    yscale = :log10,
    title = title_str,
    #colorbar_title="m^3/s",
    yticks = [1.e-1, 1.e0, 1.e1, 1.e2],
    xtickfont=20,
    ytickfont=20,
    titlefont=30,
    yguidefontsize=20,
    xguidefontsize=20,
    left_margin=10mm,
    right_margin=40mm,
    bottom_margin=10mm,
    top_margin=10mm,
    colorbar_titlefontsize=20,
    )#colorbar_title_locations=20mm)

    return h
end


function PlotDensities(output::OutputBlock,
    species_list::Vector{Species})

    n_data = output.n_data_frame

    p = plot(xlabel = "Input power [W]",
        ylabel = "n [m^-3]",
        xscale = :identity,
        #xlims = (1.e-10, time[end]),
        yscale = :log10,
        ylims = (1.e10, 1.e20)
    )

    ix = 0
    s_id = 1
    for s in species_list
        x_data = n_data[!,1]
        y_data = n_data[!,s.id + 1]
        plot!(p, x_data, y_data, w = 3, label=s.name, color=s_id,
            annotations = (x_data[end-ix], y_data[end-ix], (s.name, :bottom)) )
        ix += 7
        if ix > 32
            ix = 0
        end
        s_id += 1
    end
    plot!(p, legend=:outerleft)
    return p
end


function Plot_ChargeBalance(output::OutputBlock,
    species_list::Vector{Species})

    p = plot(xlabel = "t [ms]",
        ylabel = "n [m^-3]"
        #xscale = :log10,
        #xlims = (5.e-2, 1.e4),
        #yscale = :log10,
        #ylims = (1.e14, 1.e17)
    )

    n_data = output.n_data_frame
    time = n_data.time

    charge_balance = zeros(length(time))
    for s in species_list
        if !s.has_dens_eq
            continue
        end
        q = s.charge
        charge_dens = q * n_data[!,s.name]
        charge_balance .+= charge_dens
    end
    plot!(p, time, charge_balance, title="Charge balane", w = 2)
    return p
end


function PlotRateCoefficients(K0_df::DataFrame, reaction_list::Vector{Reaction})
    
    df_names = names(K0_df)

    c_id = 1
    p = plot()
    for r in reaction_list
        dndt = K0_df[!, r.id + 1]
        plot!(p, K0_df[!,1], dndt, label=df_names[r.id+1], w = 3, color=c_id)
        c_id += 1
    end
    plot!(p, xlim=(K0_df[1,1], K0_df[end,1]),
        xscale=:identity,
        ylim = (1.e18, 1.e22),
        yscale=:log10,
        legend=true)
    return p
end


function CompareSimulations(out1::OutputBlock, out2::OutputBlock, s_list::Vector{Species})
    # out2 is with fix O2
    n_data1 = out1.n_data_frame
    n_data2 = out2.n_data_frame

    p = plot()
    for s in s_list
        plot!(p, n_data1.time, n_data1[!,s.name],
            w = 2, linestyle=:solid, color=s.id,
            label=s.name)
        plot!(p, n_data2.time, n_data2[!,s.name],
            w = 2, linestyle=:dot, color=s.id,
            label=s.name*"(n_O2 fix)")
    end
    plot!(p, 
        xlim = (1.e-5, n_data2.time[end]),
        xscale=:log10,
        ylim = (1.e14, 1.e21),
        yscale=:log10,
        xlabel="time [s]",
        ylabel="n [m^-3]",
        legend=:topleft
    )
    return p
end


function PlotEmissionRates(K_data::DataFrame)
    # 130 nm emission from reactions
    # - "r297: O_3s -> O"
    # - "r332: e + O2 -> e + 2O"
    # - "r333: e + O2_a1Ag -> e + 2O"
    # - "r334: e + O2_b1Su -> e + 2O"
    # - "r335: e + O2_a1Ag_v -> e + 2O"
    # - "r336: e + O2_b1Su_v -> e + 2O"

    p = plot(xlabel="Total pressure [Pa]",
        ylabel="Emission rate [1/m^3/s]", xscale=:log10, xlim=(1.e-1,1e2),
        yscale=:log10, ylim=(1.e17,1.e22),
        title="Emission rates (P=1500 W)"
    )
    plot!(p, K_data[!,1], K_data."r297: O_3s -> O", label = "O_3s -> O + 130nm", w=3)
    plot!(p, K_data[!,1], K_data."r332: e + O2 -> e + 2O", label = "e + O2 -> e + 2O + 130nm", w=3)
    #plot!(p, K_data[!,1], K_data."r333: e + O2_a1Ag -> e + 2O", label = "e + O2_a1 -> e + 2O + 130nm", w=3)
    #plot!(p, K_data[!,1], K_data."r334: e + O2_b1Su -> e + 2O", label = "e + O2_b1 -> e + 2O + 130nm", w=3)
    #plot!(p, K_data[!,1], K_data."r335: e + O2_a1Ag_v -> e + 2O", label = "e + O2_a1_v -> e + 2O + 130nm", w=3)
    #plot!(p, K_data[!,1], K_data."r336: e + O2_b1Su_v -> e + 2O", label = "e + O2_b1_v -> e + 2O + 130nm", w=3)

    # Emission at 135 nm
    # - "r296: O_5s -> O"
    # - "r337: e + O2 -> e + 2O"
    # - "r338: e + O2_a1Ag -> e + 2O"
    # - "r339: e + O2_b1Su -> e + 2O"
    # - "r340: e + O2_a1Ag_v -> e + 2O"
    # - "r341: e + O2_b1Su_v -> e + 2O"
    #p = plot(xlabel="Input power [W]",
    #    ylabel="Emission rate [1/m^3/s]",
    #    yscale=:log10, ylim=(1.e16,1.e21),
    #    title="135 nm emission rates"
    #)
    plot!(p, K_data[!,1], K_data."r296: O_5s -> O", label = "O_5s -> O + 135nm", w=3)
    plot!(p, K_data[!,1], K_data."r337: e + O2 -> e + 2O", label="e + O2 -> e + 2O + 135nm", w=3)
    #plot!(p, K_data[!,1], K_data."r338: e + O2_a1Ag -> e + 2O", label="e + O2_a1 -> e + 2O + 135nm", w=3)
    #plot!(p, K_data[!,1], K_data."r339: e + O2_b1Su -> e + 2O", label="e + O2_b2 -> e + 2O + 135nm", w=3)
    #plot!(p, K_data[!,1], K_data."r340: e + O2_a1Ag_v -> e + 2O", label="e + O2_a1_v -> e + 2O + 135nm", w=3)
    #plot!(p, K_data[!,1], K_data."r341: e + O2_b1Su_v -> e + 2O", label="e + O2_b1_v -> e + 2O + 135nm", w=3)


    plot!(p, legend=:outertop)
end


function PlotEmissionRatesHeatmap(K_df::DataFrame)
    # 130 nm emission from reactions
    # - "r297: O_3s -> O"
    # - "r332: e + O2 -> e + 2O"
    # - "r333: e + O2_a1Ag -> e + 2O"
    # - "r334: e + O2_b1Su -> e + 2O"
    # - "r335: e + O2_a1Ag_v -> e + 2O"
    # - "r336: e + O2_b1Su_v -> e + 2O"
    r_list_130nm = String[]
    push!(r_list_130nm, "r297: O_3s -> O")
    push!(r_list_130nm, "r332: e + O2 -> e + 2O")
    push!(r_list_130nm, "r333: e + O2_a1Ag -> e + 2O")
    push!(r_list_130nm, "r334: e + O2_b1Su -> e + 2O")
    push!(r_list_130nm, "r335: e + O2_a1Ag_v -> e + 2O")
    push!(r_list_130nm, "r336: e + O2_b1Su_v -> e + 2O")

    # Emission at 135 nm
    # - "r296: O_5s -> O"
    # - "r337: e + O2 -> e + 2O"
    # - "r338: e + O2_a1Ag -> e + 2O"
    # - "r339: e + O2_b1Su -> e + 2O"
    # - "r340: e + O2_a1Ag_v -> e + 2O"
    # - "r341: e + O2_b1Su_v -> e + 2O"
    r_list_135nm = String[]
    push!(r_list_135nm, "r296: O_5s -> O")
    push!(r_list_135nm, "r337: e + O2 -> e + 2O")
    push!(r_list_135nm, "r338: e + O2_a1Ag -> e + 2O")
    push!(r_list_135nm, "r339: e + O2_b1Su -> e + 2O")
    push!(r_list_135nm, "r340: e + O2_a1Ag_v -> e + 2O")
    push!(r_list_135nm, "r341: e + O2_b1Su_v -> e + 2O")

    col1 = 1
    col2 = 3
    source = "/home/moe505/Documents/GM_Julia/ArO2/Param_Investigation:P_total-Power/"
    i = 1
    label_list = ["130nm", "135nm"]
    for r_list in [r_list_130nm, r_list_135nm]
        label = label_list[i]
        global row, col, matrix = Get2DMap(K_df, r_list[i], col1, col2)
        mat_tot = zeros(length(row), length(col))
        if i==1
            global mat_tot_VUV = zeros(length(row), length(col))
        end
        for r_name in r_list
            print(r_name,"\n")
            row, col, matrix = Get2DMap(K_df, r_name, col1, col2)
            mat_tot .+= matrix
            title = r_name[7:end] * " + "* label *", log10([1/m^3/s])"
            PlotHeatmap(row, col, log10.(matrix), title)
            fig_name = source * label * "_" * r_name 
            png(fig_name)
        end
        title = "Total "*label*" emission, log10([1/m^3/s])"
        PlotHeatmap(row, col, log10.(mat_tot), title)
        fig_name = source * "total_"*label*"_emission" 
        png(fig_name)

        for r_name in r_list
            row, col, matrix = Get2DMap(K_df, r_name, col1, col2)
            title = r_name[7:end] * " + "* label
            PlotHeatmap(row, col, matrix./mat_tot, title)
            fig_name = source * label * "_fraction_" * r_name 
            png(fig_name)
        end

        mat_tot_VUV .+= mat_tot
        i += 1
    end
    title = "Total VUV emission, log10([1/m^3/s])"
    PlotHeatmap(row, col, log10.(mat_tot_VUV), title)
    fig_name = source * "total_VUV_emission" 
    png(fig_name)
end

function DissociationRateHeatmap(df::DataFrame)

    col1 = 1
    col2 = 3
    source = "/home/moe505/Documents/GM_Julia/ArO2/Param_Investigation:P_total-Power/"
    row, col, matrix_O = Get2DMap(df, "O", col1, col2)
    row, col, matrix_O2 = Get2DMap(df, "O2", col1, col2)
    matrix_O_half = matrix_O .* 0.5
    dissociation = matrix_O_half ./ (matrix_O_half .+ matrix_O2)
    PlotHeatmap(row, col, dissociation, "Oxygen dissociation rate")
    fig_name = source * "oxygen_dissociation_rate" 
    png(fig_name)

end

end