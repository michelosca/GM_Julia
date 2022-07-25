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
using DataFrames, CSV
using Printf

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

        for i in range(2, length(out_list), step=1)
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

    source = "/home/moe505/Documents/GM_Julia/ArO2/ArO2_Fiebrandt_test/O_loss_test/" 
    #str_list = ["gamma_O2=0.0", "gamma_O2=0.007", "gamma_02=0.1", "gamma_02=1.0"]
    #str_list = ["gamma_02=0.007 / O_x -> O", "gamma_02=0.007 / 2O_x -> O2"]
    #str_list = ["gamma_Ar=1.0", "gamma_Ar=0.0"]
    #str_list = ["gamma_O-meta=0.0", "gamma_O-meta=0.0 / 2O_x -> O2", "gamma_O-meta=1.0", "gamma_O-meta=1.0 / 2O_x -> O2"]
    str_list = ["GM Julia","GM Julia change"]
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
            #print(r.name, "\n")
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


function PlotHeatmap(row, col, matrix, title_str, xlabel_str, ylabel_str)

    h = contourf(col,
    #h = contourf(col,
    row, matrix,
    size=(1400, 800),
    ylabel= ylabel_str,
    xlabel = xlabel_str,
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


function PlotDifferentO2fractions(main::String)
    O2_min = 0.01
    O2_max = 0.20
    step = 0.01
    n_entries = round(Int64, (O2_max - O2_min + step) / step)
    
    # 1st col: value
    # 2nd col: row (press)
    # 3rd col: col (power)
    # 4th col: O2 fraction
    abs_VUV_max = zeros(4,n_entries)  
    abs_VUV_min = zeros(4,n_entries) 
    VUV_ion_max = zeros(4,n_entries) 
    VUV_ion_min = zeros(4,n_entries) 
    VUV_O_max =   zeros(4,n_entries) 
    VUV_O_min =   zeros(4,n_entries) 

    entry = 1
    for O2_fract in range(O2_min, O2_max,step=step)
        O2_fract_str = @sprintf("%4.2f",O2_fract)
        folder = main * "O2_" * O2_fract_str * "/"
        K_file = folder * "K_vs_P_Ar_vs_P%_O2_vs_power.csv"
        param_file = folder * "param_vs_P_Ar_vs_P%_O2_vs_power.csv"
        K_dataframe = CSV.read(K_file, DataFrame)
        param_dataframe = CSV.read(param_file, DataFrame)
        max_data, min_data = PlotEmissionDataHeatmap(K_dataframe, param_dataframe, folder, 1, 3)
        # max/min data ist listed as follows:
        #  1.- 130nm emission from r297: O_3s -> O
        #  2.- 130nm emission from r332: e + O2 -> e + 2O
        #  3.- 130nm emission from r333: e + O2_a1Ag -> e + 2O 
        #  4.- 130nm emission from r334: e + O2_b1Su -> e + 2O
        #  5.- 130nm emission from r335: e + O2_a1Ag_v -> e + 2O
        #  6.- 130nm emission from r336: e + O2_b1Su_v -> e + 2O
        #  7.- 135nm emission from r296: O_5s -> O
        #  8.- 135nm emission from r337: e + O2 -> e + 2O
        #  9.- 135nm emission from r338: e + O2_a1Ag -> e + 2O
        # 10.- 135nm emission from r339: e + O2_b1Su -> e + 2O
        # 11.- 135nm emission from r340: e + O2_a1Ag_v -> e + 2O
        # 12.- 135nm emission from r341: e + O2_b1Su_v -> e + 2O
        # 13.- total emission from 130 nm
        # 14.- total emission from 135 nm
        # 15.- total VUV emission: 130 + 135 nm
        # 16.- total VUV emission / ion loss
        # 17.- total VUV emission / O loss
        abs_VUV_max[:, entry] = [max_data[15][1], max_data[15][2], max_data[15][3], O2_fract]
        abs_VUV_min[:, entry] = [min_data[15][1], min_data[15][2], min_data[15][3], O2_fract]
        VUV_ion_max[:, entry] = [max_data[16][1], max_data[16][2], max_data[16][3], O2_fract]
        VUV_ion_min[:, entry] = [min_data[16][1], min_data[16][2], min_data[16][3], O2_fract]
        VUV_O_max[:, entry]   = [max_data[17][1], max_data[17][2], max_data[17][3], O2_fract]
        VUV_O_min[:, entry]   = [min_data[17][1], min_data[17][2], min_data[17][3], O2_fract]
        val = abs_VUV_max[1,entry] 
        press = abs_VUV_max[2,entry]
        power = abs_VUV_max[3,entry] 
        @printf("O2-fract: %2.4f; Power: %15g;  Pressure: %15g, Max. abs. VUV emission %15g \n",
            O2_fract, power, press, val )
        entry += 1
    end

    #### Plots
   
    ## Max. Abs. VUV emission
    plot_data = abs_VUV_max[1,:]/1.e21
    press = abs_VUV_max[2,:]
    power = abs_VUV_max[3,:]
    O2_fract = abs_VUV_max[4,:]
    title_str = "Max. abs. VUV-emission"
    file_str = "max_VUV"
    postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)

    # Min. Abs. VUV emission
    plot_data = abs_VUV_min[1,:]/1.e21
    press = abs_VUV_min[2,:]
    power = abs_VUV_min[3,:]
    O2_fract = abs_VUV_min[4,:]
    title_str = "Min. abs. VUV-emission"
    file_str = "min_VUV"
    postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)

    # Max. VUV/ion
    plot_data = VUV_ion_max[1,:]
    press =     VUV_ion_max[2,:]
    power =     VUV_ion_max[3,:]
    O2_fract =  VUV_ion_max[4,:]
    title_str = "Max. VUV-emission/Ion-loss"
    file_str = "max_VUVion"
    postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)

    # Min. VUV/ion
    plot_data = VUV_ion_min[1,:]
    press =     VUV_ion_min[2,:]
    power =     VUV_ion_min[3,:]
    O2_fract =  VUV_ion_min[4,:]
    title_str = "Min. VUV-emission/Ion-loss"
    file_str = "min_VUVion"
    postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)

    # Max. VUV/O
    plot_data = VUV_O_max[1,:]
    press =     VUV_O_max[2,:]
    power =     VUV_O_max[3,:]
    O2_fract =  VUV_O_max[4,:]
    title_str = "Max. VUV-emission/O-loss"
    file_str = "max_VUVO"
    postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)
    
    # Min. VUV/O
    plot_data = VUV_O_min[1,:]
    press =     VUV_O_min[2,:]
    power =     VUV_O_min[3,:]
    O2_fract =  VUV_O_min[4,:]
    title_str = "Min. VUV-emission/O-loss"
    file_str = "min_VUVO"
    postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)

    return abs_VUV_max, VUV_ion_max, VUV_O_max, abs_VUV_min, VUV_ion_min, VUV_O_min
end


function postprocess_plots(press, O2_fract, power, plot_data, title_str, file_str, main)
    vuv_axis_label = "[1.e21 m^-3 s^-1]"
    press_axis_label = "P [Pa]"
    O2_axis_label = "O2 fraction [1]"
    power_axis_label = "Power [W]"
    # vs. press
    p = plot(press, plot_data, title = title_str)
    plot!(p, legend=false)
    plot!(p, xlabel=press_axis_label, ylabel=vuv_axis_label)
    png(main * file_str * "_vs_press")
    # vs. O2 fraction 
    p = plot(O2_fract, plot_data, title = title_str)
    plot!(p, xlabel=O2_axis_label, ylabel=vuv_axis_label)
    plot!(p, legend=false)
    png(main * file_str *"_vs_O2fraction")
    # vs. power 
    p = plot(power, plot_data, title = title_str)
    plot!(p, xlabel=power_axis_label, ylabel=vuv_axis_label)
    plot!(p, legend=false)
    png(main * file_str *"_vs_power")
    # 3D plots
    p = plot3d( press, O2_fract, plot_data, title=title_str)
    plot3d!(p, legend=false)
    plot3d!(p, xlabel=press_axis_label, ylabel=O2_axis_label, zlabel=vuv_axis_label)
    png(main * file_str *"_vs_press_vs_O2fraction")

    # Generate a heatmap plot
    press_min = 0.1 
    press_max = 100.0
    press_steps = 95
    press_arr = [10^y for y in range(log10(press_min), log10(press_max), length=press_steps)]
    O2f_min = 0.01
    O2f_max = 0.2
    O2f_steps = 20
    O2f_arr = range(O2f_min,O2f_max, length=O2f_steps)
    data_mat = zeros(press_steps, O2f_steps)

    len_data = length(press)
    for i in range(1, length=len_data)
        p = press[i]
        O2f = O2_fract[i]
        data_point = plot_data[i]
        mat_p_ix = 0
        mat_O2f_ix = 0
        # Check position of pressure in matrix
        for p_ix in range(1,length=press_steps)
            p_val = press_arr[p_ix]
            p_half_step_min = 0.0
            p_half_step_max = 0.0
            if p_ix == 1
                p_half_step_max = (press_arr[p_ix+1]-p_val)*0.5
                p_half_step_min = p_half_step_max
            elseif p_ix == press_steps
                p_half_step_min = (p_val-press_arr[p_ix-1])*0.5
                p_half_step_max = p_half_step_min
            else
                p_half_step_min = (p_val - press_arr[p_ix-1])*0.5
                p_half_step_max = (press_arr[p_ix+1]-p_val)*0.5
            end
            p_min = p_val - p_half_step_min
            p_max = p_val + p_half_step_max
            if p <= p_max && p >= p_min
                mat_p_ix = p_ix
                break
            end
        end
        # Check position of O2 fraction in matrix
        for O2f_ix in range(1,length=O2f_steps)
            O2f_val = O2f_arr[O2f_ix]
            O2f_half_step = 0.01 * 0.5
            O2f_min = O2f_val - O2f_half_step
            O2f_max = O2f_val + O2f_half_step
            #print("Value tested ", O2f, "boundaries ", O2f_min, " ", O2f_max,"\n")
            if O2f <= O2f_max && O2f >= O2f_min
                #print("VAlue found!\n")
                mat_O2f_ix = O2f_ix
                break
            end
        end
        data_mat[mat_p_ix, mat_O2f_ix] = data_point
    end
    heatmap(O2f_arr, press_arr, data_mat,
    size=(1400, 800),
    ylabel= "Pressure [W]",     # row data label
    xlabel = "O2 fraction [1]",  # row data label
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
    png(main * file_str *"_vs_press_vs_O2fraction_heatmap")

    open(main*file_str*".dat","w") do file
        @printf(file, "%25s %25s %25s %25s\n","Pressure [Pa]", "Power [W]", "O2-fraction [1]", title_str)
        for i in range(1,length=length(plot_data))
            @printf(file, "%25g %25g %25g %25g\n", press[i], power[i], O2_fract[i], plot_data[i])
        end
    end
end


function PlotEmissionDataHeatmap(K_df::DataFrame, param_df::DataFrame,
    source::String, col1::Int64, col2::Int64)

    xlabel = "Power [W]"
    ylabel = "Pressure [Pa]"
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

    # Generate data for O-atom loss rate
    row, col, Oflux_mat = Get2DMap(K_df, "r281: 2O -> O2", col1, col2)

    # Generate matrix for ion-flux loss rate
    row, col, nOion_mat  = Get2DMap(param_df, "O+_flux",  col1, col2)
    row, col, nO2ion_mat = Get2DMap(param_df, "O2+_flux", col1, col2)
    row, col, nO3ion_mat = Get2DMap(param_df, "O3+_flux", col1, col2)
    row, col, nO4ion_mat = Get2DMap(param_df, "O4+_flux", col1, col2)
    row, col, nArion_mat = Get2DMap(param_df, "Ar+_flux", col1, col2)
    ionflux_mat = nOion_mat .+ nO2ion_mat .+ nO3ion_mat .+ nO4ion_mat .+ nArion_mat
    r = 0.2
    L = 0.2
    area = pi * r^2 * 2 + 2.0*pi*r * L
    volume = pi * r^2 * L
    ionflux_mat *= area / volume

    max_data = Tuple{Float64, Float64, Float64}[]
    min_data = Tuple{Float64, Float64, Float64}[]
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
            #print(r_name,"\n")
            row, col, matrix = Get2DMap(K_df, r_name, col1, col2)

            # Get max/min data and store it in list
            val, row_ix, col_ix = GetMaxValue(matrix)
            push!(max_data, (val, row[row_ix], col[col_ix]))
            val, row_ix, col_ix = GetMinValue(matrix)
            push!(min_data, (val, row[row_ix], col[col_ix]))

            # Add matrix to total emission matrix
            mat_tot .+= matrix

            # Plot and save figure
            title = r_name[7:end] * " + "* label *", log10([1/m^3/s])"
            PlotHeatmap(row, col, log10.(matrix), title, xlabel, ylabel)
            fig_name = source * label * "_" * r_name[1:4] 
            fig_name = join(map(x -> isspace(fig_name[x]) ? "" : fig_name[x], 1:length(fig_name)))
            png(fig_name)
        end
        # Get max/min data and store it in list
        val, row_ix, col_ix = GetMaxValue(mat_tot)
        push!(max_data, (val, row[row_ix], col[col_ix]))
        val, row_ix, col_ix = GetMinValue(mat_tot)
        push!(min_data, (val, row[row_ix], col[col_ix]))
        
        # Plot total emission and save figure
        title = "Total "*label*" emission, log10([1/m^3/s])"
        PlotHeatmap(row, col, log10.(mat_tot), title, xlabel, ylabel)
        fig_name = source * "total_"*label*"_emission" 
        fig_name = join(map(x -> isspace(fig_name[x]) ? "" : fig_name[x], 1:length(fig_name)))
        png(fig_name)

        # Generate fraction plots and save figures
        for r_name in r_list
            row, col, matrix = Get2DMap(K_df, r_name, col1, col2)
            title = r_name[7:end] * " + "* label
            PlotHeatmap(row, col, matrix./mat_tot, title, xlabel, ylabel)
            fig_name = source * label * "_fraction_" * r_name[1:4] 
            fig_name = join(map(x -> isspace(fig_name[x]) ? "" : fig_name[x], 1:length(fig_name)))
            png(fig_name)
        end

        mat_tot_VUV .+= mat_tot
        i += 1
    end

    ### TOTAL VUV EMISSION ####################################################
    # Get max/min data and store it in list
    val, row_ix, col_ix = GetMaxValue(mat_tot_VUV)
    push!(max_data, (val, row[row_ix], col[col_ix]))
    val, row_ix, col_ix = GetMinValue(mat_tot_VUV)
    push!(min_data, (val, row[row_ix], col[col_ix]))

    # Generate total VUV emission and save plot
    title = "Total VUV emission, log10([1/m^3/s])"
    PlotHeatmap(row, col, log10.(mat_tot_VUV), title, xlabel, ylabel)
    fig_name = source * "total_VUV_emission" 
    png(fig_name)

    ### VUV-EMISSION/ION-LOSS RATE #############################################
    # Get max/min data and store it in list
    vuv_ion_rate = mat_tot_VUV ./ ionflux_mat
    val, row_ix, col_ix = GetMaxValue(vuv_ion_rate)
    push!(max_data, (val, row[row_ix], col[col_ix]))
    val, row_ix, col_ix = GetMinValue(vuv_ion_rate)
    push!(min_data, (val, row[row_ix], col[col_ix]))

    # Generate total VUV emission and save plot
    ### Log10 scale
    title = "VUV emission / ion-loss rate [log10()]"
    PlotHeatmap(row, col, log10.(vuv_ion_rate), title, xlabel, ylabel)
    fig_name = source * "VUV_ionloss_rate_log10" 
    png(fig_name)
    ### Linear scale
    title = "VUV emission / ion-loss rate"
    PlotHeatmap(row, col, vuv_ion_rate, title, xlabel, ylabel)
    fig_name = source * "VUV_ionloss_rate" 
    png(fig_name)
    
    ### VUV-EMISSION/O-LOSS RATE #############################################
    # Get max/min data and store it in list
    vuv_O_rate = mat_tot_VUV ./ Oflux_mat
    val, row_ix, col_ix = GetMaxValue(vuv_O_rate)
    push!(max_data, (val, row[row_ix], col[col_ix]))
    val, row_ix, col_ix = GetMinValue(vuv_O_rate)
    push!(min_data, (val, row[row_ix], col[col_ix]))

    # Generate total VUV emission and save plot
    ### Linear scale
    title = "VUV emission / O-loss rate"
    PlotHeatmap(row, col, vuv_O_rate, title, xlabel, ylabel)
    fig_name = source * "vuv_Oloss_rate" 
    png(fig_name)
    ### Log10 scale
    title = "VUV emission / O-loss rate [log10()]"
    PlotHeatmap(row, col, log10.(vuv_O_rate), title, xlabel, ylabel)
    fig_name = source * "vuv_Oloss_rate_log10" 
    png(fig_name)


    return max_data, min_data
end


function DissociationRateHeatmap(source::String, col1::Int64, col2::Int64)

    n_df = CSV.read(source * "n_vs_P_Ar_vs_P%_O2_vs_power.csv", DataFrame)
    row, col, matrix_O = Get2DMap(n_df, "O", col1, col2)
    row, col, matrix_O2 = Get2DMap(n_df, "O2", col1, col2)
    matrix_O_half = matrix_O .* 0.5
    dissociation = matrix_O_half ./ (matrix_O_half .+ matrix_O2)
    title_str = "Oxygen dissociation rate"
    xlabel_str = "Power [W]"
    ylabel_str = "Pressure [Pa]"
    h = contourf(col, row, dissociation,
        ylabel= ylabel_str, xlabel = xlabel_str, title = title_str,
        #colorbar_title="m^3/s",
        size=(1400, 800),
        xscale = :identity, yscale = :log10,
        yticks = [1.e-1, 1.e0, 1.e1, 1.e2],
        xtickfont=20, ytickfont=20, titlefont=30,
        yguidefontsize=20, xguidefontsize=20, colorbar_titlefontsize=20,
        left_margin=10mm, right_margin=40mm, bottom_margin=10mm, top_margin=10mm
    ) #colorbar_title_locations=20mm)
    fig_name = source * "oxygen_dissociation_rate" 
    png(fig_name)

end


function DensityHeatmap(source::String, col1::Int64, col2::Int64, species::String)

    n_df = CSV.read(source * "n_vs_P_Ar_vs_P%_O2_vs_power.csv", DataFrame)
    row, col, matrix = Get2DMap(n_df, species, col1, col2)
    title_str = species * " number density [log10(m^-3)]" 
    xlabel_str = "Power [W]"
    ylabel_str = "Pressure [Pa]"
    h = contourf(col, row, log10.(matrix),
        ylabel= ylabel_str, xlabel = xlabel_str, title = title_str,
        #colorbar_title="m^3/s",
        size=(1400, 800),
        xscale = :identity, yscale = :log10,
        yticks = [1.e-1, 1.e0, 1.e1, 1.e2],
        xtickfont=20, ytickfont=20, titlefont=30,
        yguidefontsize=20, xguidefontsize=20, colorbar_titlefontsize=20,
        left_margin=10mm, right_margin=40mm, bottom_margin=10mm, top_margin=10mm
    ) #colorbar_title_locations=20mm)
    fig_name = source * species * "_number_density" 
    png(fig_name)

end


function TempHeatmap(source::String, col1::Int64, col2::Int64)

    T_df = CSV.read(source * "T_vs_P_Ar_vs_P%_O2_vs_power.csv", DataFrame)
    row, col, matrix = Get2DMap(T_df, "e", col1, col2)
    title_str = "electron temperature [eV]" 
    xlabel_str = "Power [W]"
    ylabel_str = "Pressure [Pa]"
    h = contourf(col, row, matrix * K_to_eV,
        ylabel= ylabel_str, xlabel = xlabel_str, title = title_str,
        #colorbar_title="m^3/s",
        size=(1400, 800),
        xscale = :identity, yscale = :log10,
        yticks = [1.e-1, 1.e0, 1.e1, 1.e2],
        xtickfont=20, ytickfont=20, titlefont=30,
        yguidefontsize=20, xguidefontsize=20, colorbar_titlefontsize=20,
        left_margin=10mm, right_margin=40mm, bottom_margin=10mm, top_margin=10mm
    ) #colorbar_title_locations=20mm)
    fig_name = source * "e_temperature" 
    png(fig_name)

end


function GetMaxValue(matrix::Array{Float64,2})
    max_loc = findall( x -> x == maximum(matrix), matrix)
    row = max_loc[1][1]
    col = max_loc[1][2]
    max_val = matrix[row, col]
    return max_val, row, col
end


function GetMinValue(matrix::Array{Float64,2})
    max_loc = findall( x -> x == minimum(matrix), matrix)
    row = max_loc[1][1]
    col = max_loc[1][2]
    min_val = matrix[row, col]
    return min_val, row, col
end

function FindReactionPathways(main_species::Species,
    reaction_list::Vector{Reaction}, species_list::Vector{Species},
    sID::SpeciesID, source::String)

    # Load rate coefficient data frame
    K_df = CSV.read(source * "K_vs_P_Ar_vs_P%_O2_vs_power.csv", DataFrame)
    ## Generate dummy data in order to get number of rows and columns
    row, col, matrix = Get2DMap(K_df, "r1: e + O -> 2e + O+", 1, 3)
    nrows = length(row)
    ncols = length(col)

    # Make a directory where the figures are going to be saved
    fold_container = main_species.name * "_pathways"
    try
        mkdir(source * fold_container)
    catch
        print("*** ERROR *** production_rates folder already exists\n")
        return
    end

    # Set chain counter
    pathway_chain = 0
    
    # Start loop: find reaction pathways for the 'main_species'
    #  - find reactions productin 'main_species'
    #  - selects source species
    #  - tracks these new species and the reactions that generate them
    loop_flag = true
    species_id_list = Int64[main_species.id]
    species_check_list = species_id_list
    while loop_flag
        loop_flag = false

        # New species found are added to the buffer list and added later to 'species_id_list'
        species_id_buffer = Int64[]

        # Loop over the species list
        for s in species_id_list

            global species_fold = @sprintf("%2.2i_%s", pathway_chain, species_list[s].name)
            mkdir(source * fold_container * "/" * species_fold)

            # Find reactions that produce the species
            reaction_products = FindProductSpecies(reaction_list, species_list[s])

            # Get total reaction rate producing species 's'
            K_total = zeros(nrows, ncols)
            for r in reaction_products
                K_name = "r" * string(r.id) * ": "*r.name
                row, col, K_data = Get2DMap(K_df, K_name, 1, 3)
                K_total .+= K_data
            end

            # Plot reaction data that produces species 's' 
            for r in reaction_products
                # Load data into matrix
                K_name = "r" * string(r.id) * ": "*r.name
                row, col, K_data = Get2DMap(K_df, K_name, 1, 3)

                # Reaction rate over total rate coefficient
                K_rate = K_data ./ K_total

                # Data is only relevant if it is more than 10% of the total rate coefficient
                plot_flag = maximum(K_rate) >= 0.10
                if plot_flag
                    # Save product species of the reaction. Only if:
                    # - Is not an electron species
                    # - It has not already in the species_id_list
                    for s_react_id in r.reactant_species
                        # Skip electron species
                        if s_react_id == sID.electron
                            continue
                        end

                        # Check wether the species is already in the list
                        list_flag = ( findall(x-> x == s_react_id, species_check_list)==Int64[] ) 
                        buff_flag = ( findall(x-> x == s_react_id, species_id_buffer)==Int64[] ) 
                        push_flag = list_flag && buff_flag

                        # Add species id only if not found in the species-list yet
                        if push_flag
                            push!(species_id_buffer, s_react_id)
                        end
                    end # end loop over involved species in reaction

                    # Plot data
                    title_str = r.name * " : " * species_list[s].name * " production rate"
                    xlabel_str = "Power [W]"
                    ylabel_str = "Pressure [Pa]"
                    h = contourf(col, row, K_rate,
                        ylabel= ylabel_str, xlabel = xlabel_str, title = title_str,
                        #colorbar_title="m^3/s",
                        size=(1400, 800),
                        xscale = :identity, yscale = :log10,
                        yticks = [1.e-1, 1.e0, 1.e1, 1.e2],
                        xtickfont=20, ytickfont=20, titlefont=30,
                        yguidefontsize=20, xguidefontsize=20, colorbar_titlefontsize=20,
                        left_margin=10mm, right_margin=40mm, bottom_margin=10mm, top_margin=10mm
                    ) #colorbar_title_locations=20mm)
                    file_name = species_list[s].name * "_r" * string(r.id) * "_" * "production_rate"
                    png(source * fold_container * "/" * species_fold * "/" * file_name)
                end

            end # end loop reaction_products list 
        end # end loop over species_id_list

        # Update species_id_list with species found
        species_check_list = [species_check_list;species_id_buffer]
        species_id_list = species_id_buffer
        if species_id_buffer != Int64[]
            loop_flag = true
        end
        pathway_chain += 1
    end
end

end