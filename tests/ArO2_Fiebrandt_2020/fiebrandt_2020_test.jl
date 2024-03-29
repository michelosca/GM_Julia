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

module fiebrandt_2020_test 

using Plots
using Plots.PlotMeasures
using LaTeXStrings
using SharedData: Species, OutputBlock, Reaction, SpeciesID, System
using SharedData: K_to_eV, kb
using DataFrames, CSV
using Printf

global ix_ref = 299 # reaction reference is first VUV emitting e + O -> e + O_5s

function CRM_test(n_data::DataFrame, T_data::DataFrame, K_data::DataFrame, output_fold::String)

    # Fiebrandt's data
    O2_fract_fieb20 = [2.5/102.5, 5/105, 7.5/107.5,10/110, 20/120, 1]
    T_g_fieb2020_K = [622,683,717,739,787,700]

    # GM data
    O2_fract_GM = n_data."P%_O2"

    # Plot params
    plot_font = "Computer modern"
    default(fontfamily=plot_font
        , linewidth=2
        , framestyle=:box
        , shape = :none
    )
    GM_label = "GM-Julia"
    F20_label = "Fiebrandt 2020"
    F17_label = "Fiebrandt 2017"

    # Electron temeperature 
    Te_fieb2020_eV = [1.9, 2.2, 2.5, 2.6, 2.7, 2]
    Te_eV = T_data.e * K_to_eV
    p = plot(O2_fract_GM, Te_eV
        , color = :black
        , label = GM_label
        , legend = :bottomright
        , ylims = (0, 4)
        , ylabel = L"\textrm{\mathrm{T_e}\ /\ eV}"
    )
    scatter!(p, O2_fract_fieb20, Te_fieb2020_eV
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "Te.png")

    # Electron density
    ne_fieb2020_e16 = [52.10, 33.00, 21.23, 16.29, 10.84, 5.18]
    ne_e16 = n_data.e * 1.e-16
    p = plot(O2_fract_GM, ne_e16
        , color = :black
        , label = GM_label
        , legend = :topright
        , yscale=:log10
        , ylims = (1, 1.e3)
        , ylabel = L"\textrm{\mathrm{n_e}\ /\ \mathrm{m^{_3}10^{16}}}"
    )
    scatter!(p, O2_fract_fieb20, ne_fieb2020_e16
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "ne.png")


    # Ar_m densities
    Arm_fieb2020_e16 = [2.2+0.23, 2.8+0.28, 2.4+0.16, 1.9+0.08, 0.31+0.01, 0.0]
    Arm_fieb2017_e16 = [3.07235, 3.48848, 2.98848, 2.44674, 1.853219, 0.0]
    Arm_e16 = n_data.Ar_m * 1.e-16
    p = plot(O2_fract_GM, Arm_e16
        , color = :black
        , label = GM_label
        , legend = :bottomleft
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{Ar^m}}\ /\ \mathrm{m^{_3}10^{16}}}"
    )
    scatter!(p, O2_fract_fieb20, Arm_fieb2020_e16
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    scatter!(p, O2_fract_fieb20, Arm_fieb2017_e16
        , shape = :utriangle
        , markercolor = :white
        , markerstrokecolor = :red
        , markersize = 6
        , label = F17_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "Ar_m.png")

    # Ar_m densities
    Arr_fieb2020_e16 = [0.68+0.55, 0.82+0.67, 0.61+0.47, 0.46+0.32, 0.07+0.05, 0.0]
    Arr_fieb2017_e16 = [ 1.59728, 1.68462, 1.29639, 1.00538, 0.64125, 0.0]
    Arr_e16 = n_data.Ar_r * 1.e-16
    p = plot(O2_fract_GM, Arr_e16
        , color = :black
        , label = GM_label
        , legend = :topright
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{Ar^r}}\ /\ \mathrm{m^{_3} 10^{^16}}}"
    )
    scatter!(p, O2_fract_fieb20, Arr_fieb2020_e16
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    scatter!(p, O2_fract_fieb20, Arr_fieb2017_e16
        , shape = :utriangle
        , markercolor = :white
        , markerstrokecolor = :red
        , markersize = 6
        , label = F17_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "Ar_r.png")

    # VUV emission rates
    I_130_fieb2020_e18 = [450, 510, 460, 400, 260,82]
    I_135_fieb2020_e18 = [29,37,40,35,18,1.3]
    I_130_gm_e18 = Get130nmEmission(K_data) * 1.e-18
    I_135_gm_e18 = Get135nmEmission(K_data) * 1.e-18
    p = plot(O2_fract_GM, I_130_gm_e18 
        , color = :black
        , label = GM_label * " 130 nm" 
        , legend = :bottomleft
        , yscale=:log10
        , ylims = (1.e0, 1.e3)
        , ylabel = L"\textrm{\mathrm{I}\ /\ \mathrm{m^{_3}s^{-1} 10^{18}}}"
    )
    plot!(p, O2_fract_GM, I_135_gm_e18 
        , ls = :dot
        , color = :black
        , label = GM_label * " 135 nm" 
    )
    scatter!(p, O2_fract_fieb20, I_130_fieb2020_e18
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label * " 130 nm"
    )
    scatter!(p, O2_fract_fieb20, I_135_fieb2020_e18
        , shape = :square
        , markercolor = :white
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label * " 135 nm"
    )
    CustomizePlot!(p)
    savefig(output_fold * "VUV.png")

    # 777nm emission line
    I_777_fieb2020_e18 = [520,570,470,340,150,32]
    I_777_gm_e18 = Get777nmEmission(K_data) * 1.e-18
    p = plot(O2_fract_GM, I_777_gm_e18 
        , color = :black
        , label = GM_label
        , legend = :bottomleft
        , yscale=:log10
        , ylims = (1.e1, 1.e3)
        , ylabel = L"\textrm{\mathrm{I_{777}}\ /\ \mathrm{m^{_3}s^{-1} 10^{18}}}"
    )
    scatter!(p, O2_fract_fieb20, I_777_fieb2020_e18
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "I777.png")

    # SELF-ABSORPTION
    # 777nm self-absorption coefficient
    A_777nm = 3.69e7
    n_O5p = n_data.O_5p
    gamma_777nm_1 = K_data[!,@sprintf("r%i: O_5p -> O_5s", ix_ref+258)]/0.47 ./n_O5p/A_777nm
    gamma_777nm_2 = K_data[!,@sprintf("r%i: O_5p -> O_5s", ix_ref+259)]/0.33 ./n_O5p/A_777nm
    gamma_777nm_3 = K_data[!,@sprintf("r%i: O_5p -> O_5s", ix_ref+260)]/0.2  ./n_O5p/A_777nm
    gamma_777nm = (gamma_777nm_1  + gamma_777nm_2 + gamma_777nm_3) ./ n_O5p / A_777nm
    p = plot(O2_fract_GM, gamma_777nm_1
        , color = :black
        , label = GM_label * "-1"
        , ylim = (0,1.1)
        , ls = :solid
    )
    plot!(p, O2_fract_GM, gamma_777nm_2
        , color = :black
        , label = GM_label * "-2"
        , ls = :dash
    )
    plot!(p, O2_fract_GM, gamma_777nm_3
        , color = :black
        , label = GM_label * "-3"
        , ls = :dot
    )
    SA_777_fieb2020 = [0.47, 0.43, 0.42, 0.45, 0.63, 0.95]
    scatter!(p, O2_fract_fieb20, SA_777_fieb2020 
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "selfabsorption_I777.png")

    # 130nm self-absorption coefficient
    A_130nm_1 = 3.41e8
    A_130nm_2 = 2.03e8
    A_130nm_3 = 6.76e7
    n_O3s = n_data.O_3s
    gamma_130nm_1 = K_data[!, @sprintf("r%i: O_3s -> O", ix_ref+255)]./n_O3s/A_130nm_1*3
    gamma_130nm_2 = K_data[!, @sprintf("r%i: O_3s -> O", ix_ref+256)]./n_O3s/A_130nm_2*3
    gamma_130nm_3 = K_data[!, @sprintf("r%i: O_3s -> O", ix_ref+257)]./n_O3s/A_130nm_3*3
    gamma_130nm = (gamma_130nm_1  + gamma_130nm_2 + gamma_130nm_3) 
    p = plot(O2_fract_GM, gamma_130nm_1
        , color = :black
        , label = GM_label * "-1"
        , ylim = (0,1.1)
        , ls = :solid
    )
    plot!(O2_fract_GM, gamma_130nm_2
        , color = :black
        , label = GM_label * "-2"
        , ls = :dot
    )
    plot!(O2_fract_GM, gamma_130nm_3
        , color = :black
        , label = GM_label * "-3"
        , ls = :dash
    )
    SA_130_fieb2020 = [0.125, 0.06, 0.05, 0.04, 0.03, 0.02]
    scatter!(p, O2_fract_fieb20, SA_130_fieb2020 
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "selfabsorption_I130.png")

    # 844nm self-absorption
    A_844nm =3.22e7
    n_O3p = n_data.O_3p
    gamma_844nm_1 = K_data[!,@sprintf("r%i: O_3p -> O_3s", ix_ref+261)]/0.11./n_O3p/A_844nm
    gamma_844nm_2 = K_data[!,@sprintf("r%i: O_3p -> O_3s", ix_ref+262)]/0.56./n_O3p/A_844nm
    gamma_844nm_3 = K_data[!,@sprintf("r%i: O_3p -> O_3s", ix_ref+263)]/0.33./n_O3p/A_844nm
    gamma_844nm = (gamma_844nm_1 + gamma_844nm_2 + gamma_844nm_3)
    p = plot(O2_fract_GM, gamma_844nm_1
        , color = :black
        , label = GM_label * "-1"
        , ylim = (0,1.1)
        , ls = :solid
        , legend = :bottomleft
    )
    plot!(O2_fract_GM, gamma_844nm_2
        , color = :black
        , label = GM_label * "-2"
        , ls = :dot
    )
    plot!(O2_fract_GM, gamma_844nm_3
        , color = :black
        , label = GM_label * "-3"
        , ls = :dash
    )
    SA_844_fieb2020 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    scatter!(p, O2_fract_fieb20, SA_844_fieb2020 
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "selfabsorption_I844.png")


    # Dissociation rate
    dissociation_fieb2020_percent = [29.5, 30.8, 28.3, 25.1, 17.1, 3.1]
    nO_total = n_data.O + n_data.O_1d + n_data.O_1s + n_data.O_3s + n_data.O_5s + n_data.O_3p + n_data.O_5p
    nO2 = n_data.O2 .+ n_data.O2_a1Ag .+ n_data.O2_b1Su #.+ n_data.O2_a1Ag_v .+ n_data.O2_b1Su_v
    diss = nO_total * 0.5 ./ (0.5 * nO_total .+ nO2) * 100
    p = plot(O2_fract_GM, diss 
        , color = :black
        , label = GM_label
        , legend = :bottomleft
        , ylims = (0,40)
        , ylabel = L"\textrm{\textrm{diss.\ }\mathrm{O_2}\ /\ \%}"
    )
    scatter!(p, O2_fract_fieb20, dissociation_fieb2020_percent
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "dissO2.png")

    # O2 Molecules oxygen
    nO2_fieb2020_e20 = [0.1, 0.17, 0.25, 0.33, 0.64, 5.02]
    nO2_e20 = nO2 * 1.e-20
    p = plot(O2_fract_GM, nO2_e20
        , color = :black
        , label = GM_label
        , legend = :topleft
        , yscale=:log10
        , ylims = (6.e-2, 1.e0)
        , ylabel = L"\textrm{\mathrm{n_{O_2}}\ /\ \mathrm{m^{_3}10^{20}}}"
    )
    scatter!(p, O2_fract_fieb20, nO2_fieb2020_e20
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO2.png")


    # Atomic oxygen: GROUND STATE
    nO_GM_e19 = n_data.O * 1.e-19
    nO_fieb2020_e19 = [0.41, 0.93, 1.4, 1.7, 2.3, 3.15]
    p = plot(O2_fract_GM, nO_GM_e19
        , color = :black
        , label = GM_label
        , legend = :topleft
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{O}}\ /\ \mathrm{m^{_3}10^{19}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_fieb2020_e19
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_ground.png")

    # Atomic oxygen: 1D
    nO1d_GM_e18 = n_data.O_1d * 1.e-18
    nO_1d_fieb2020_e18 = [4.1, 6.1, 5.9, 5.2, 3.1, 0.23]
    p = plot(O2_fract_GM, nO1d_GM_e18
        , color = :black
        , label = GM_label
        , legend = :bottomright
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{O(^1D)}}\ /\ \mathrm{m^{_3}10^{18}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_1d_fieb2020_e18
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_1D.png")

    # Atomic oxygen: 1S
    nO1s_GM_e17 = n_data.O_1s * 1.e-17
    nO_1s_fieb2020_e17 = [1.3, 1.5, 0.79, 0.73, 0.54,0.11]
    p = plot(O2_fract_GM, nO1s_GM_e17
        , color = :black
        , label = GM_label
        , legend = :topright
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{O(^1S)}}\ /\ \mathrm{m^{_3}10^{18}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_1s_fieb2020_e17
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_1S.png")

    # Atomic oxygen: 3S
    nO3s_GM_e13 = n_data.O_3s * 1.e-13
    nO_3s_fieb2020_e13 = [0.6, 1.4, 1.7, 1.8, 1.5, 0.66]
    p = plot(O2_fract_GM, nO3s_GM_e13
        , color = :black
        , label = GM_label
        , legend = :topright
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{O(^3S)}}\ /\ \mathrm{m^{_3}10^{13}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_3s_fieb2020_e13
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_3S.png")
    
    # Atomic oxygen: 5S
    nO5s_GM_e15 = n_data.O_5s * 1.e-15
    nO_5s_fieb2020_e15 = [5.2, 6.7, 7.2, 6.3, 3.2, 0.23]
    p = plot(O2_fract_GM, nO5s_GM_e15
        , color = :black
        , label = GM_label
        , legend = :bottomright
        , yscale=:log10
        , ylims = (1.e0, 1.e2)
        , ylabel = L"\textrm{\mathrm{n_{O(^5S)}}\ /\ \mathrm{m^{_3}10^{13}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_5s_fieb2020_e15
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_5S.png")

    # Atomic oxygen: 3P
    nO3p_GM_e12 = n_data.O_3p * 1.e-12
    nO_3p_fieb2020_e12 = [4.7, 6.9, 9.0, 9.2, 7.2, 1.1]
    p = plot(O2_fract_GM, nO3p_GM_e12
        , color = :black
        , label = GM_label
        , legend = :bottomright
        , yscale=:log10
        , ylims = (1.e0, 1.e2)
        , ylabel = L"\textrm{\mathrm{n_{O(^3P)}}\ /\ \mathrm{m^{_3}10^{13}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_3p_fieb2020_e12
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_3P.png")

    # Atomic oxygen: 5P
    nO5p_GM_e13 = n_data.O_5p * 1.e-13
    nO_5p_fieb2020_e13 = [2.9, 3.5, 3.0, 2.0, 0.65, 0.09]
    p = plot(O2_fract_GM, nO5p_GM_e13
        , color = :black
        , label = GM_label
        , legend = :topright
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{O(^5P)}}\ /\ \mathrm{m^{_3}10^{13}}}"
    )
    scatter!(p, O2_fract_fieb20, nO_5p_fieb2020_e13
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_5P.png")

    # Atomic oxygen: GROUND STATE
    nOall_fieb2020_e19 = [0.83, 1.55, 2.0, 2.23, 2.62, 3.17]
    nOall_GM_e19 = (n_data.O + n_data.O_1d + n_data.O_1s + n_data.O_3s + n_data.O_5s + n_data.O_3p + n_data.O_5p) * 1.e-19
    p = plot(O2_fract_GM, nOall_GM_e19
        , color = :black
        , label = GM_label
        , legend = :topleft
        , yscale=:log10
        , ylims = (1.e-1, 1.e1)
        , ylabel = L"\textrm{\mathrm{n_{O}}\ /\ \mathrm{m^{_3}10^{19}}}"
    )
    scatter!(p, O2_fract_fieb20, nOall_fieb2020_e19
        , shape = :cross
        , markercolor = :red
        , markerstrokecolor = :red
        , markersize = 6
        , label = F20_label
    )
    CustomizePlot!(p)
    savefig(output_fold * "nO_all.png")

    return p
end

function CustomizePlot!(p)
    o2_min = 0
    o2_max = 0.21
    plot!(p
        , frame=:box
        , size = (500*2.5,300*2.5)
        , tickfontsize=15
        , guidefontsize=15
        , legendfontsize=10
        , legendtitlefontsize=10
        , xlabel = L"\textrm{\mathrm{O_2}\ fraction}"
        , xlims = (o2_min, o2_max)
        , xticks = ([0,0.05, 0.1, 0.15, 0.2],["0","0.05","0.1","0.15","0.2"])
        , bottommargin = 3mm
        , leftmargin = 8mm
        , topmargin = 2mm
    )
end

function Get777nmEmission(K_data::DataFrame)
    
    I_777nm =  K_data[!,@sprintf("r%i: O_5p -> O_5s", ix_ref+258)]
    I_777nm += K_data[!,@sprintf("r%i: O_5p -> O_5s", ix_ref+259)]
    I_777nm += K_data[!,@sprintf("r%i: O_5p -> O_5s", ix_ref+260)]
    I_777nm += K_data[!,@sprintf("r%i: e + O -> e + O_5s", ix_ref)]
    I_777nm += K_data[!,@sprintf("r%i: e + O2 -> e + O + O_5s", ix_ref+7)]
    I_777nm += K_data[!,@sprintf("r%i: e + O2_a1Ag -> e + O + O_5s", ix_ref+8)]
    I_777nm += K_data[!,@sprintf("r%i: e + O2_b1Su -> e + O + O_5s", ix_ref+9)]
    #I_777nm += K_data[!,@sprintf("r%i: e + O2_a1Ag_v -> e + O + O_5s", ix_ref+13)]
    #I_777nm += K_data[!,@sprintf("r%i: e + O2_b1Su_v -> e + O + O_5s", ix_ref+14)]
    return I_777nm
end

function Get844nmEmission(K_data::DataFrame, all::Bool)
    
    I_844nm =  K_data[!,@sprintf("r%i: O_3p -> O_3s", ix_ref+261)]
    I_844nm += K_data[!,@sprintf("r%i: O_3p -> O_3s", ix_ref+262)]
    I_844nm += K_data[!,@sprintf("r%i: O_3p -> O_3s", ix_ref+263)]
    I_844nm += K_data[!,@sprintf("r%i: e + O2 -> e + O + O_3s", ix_ref+10)]
    I_844nm += K_data[!,@sprintf("r%i: e + O2_a1Ag -> e + O + O_3s", ix_ref+11)]
    I_844nm += K_data[!,@sprintf("r%i: e + O2_b1Su -> e + O + O_3s", ix_ref+12)]
    #I_844nm += K_data[!,@sprintf("r%i: e + O2_a1Ag_v -> e + O + O_3s", ix_ref+18)]
    #I_844nm += K_data[!,@sprintf("r%i: e + O2_b1Su_v -> e + O + O_3s", ix_ref+19)]

    return I_844nm
end

function Get577nmEmission(K_data::DataFrame)
    
    I_544nm = K_data[!, @sprintf("r%i: O_1s -> O_1d",ix_ref+252)]
    return I_544nm
end

function Get130nmEmission(K_data::DataFrame)
    
    I_130nm = K_data[!, @sprintf("r%i: O_3s -> O", ix_ref+255)]
    I_130nm += K_data[!, @sprintf("r%i: O_3s -> O", ix_ref+256)]
    I_130nm += K_data[!, @sprintf("r%i: O_3s -> O", ix_ref+257)]
    I_130nm += K_data[!,@sprintf("r%i: e + O2 -> e + 2O", ix_ref+1)]
    I_130nm += K_data[!,@sprintf("r%i: e + O2_a1Ag -> e + 2O", ix_ref+2)]
    I_130nm += K_data[!,@sprintf("r%i: e + O2_b1Su -> e + 2O", ix_ref+3)]
    #I_130nm += K_data[!,@sprintf("r%i: e + O2_a1Ag_v -> e + 2O", ix_ref+4)]
    #I_130nm += K_data[!,@sprintf("r%i: e + O2_b1Su_v -> e + 2O", ix_ref+5)]
    return I_130nm
end

function Get135nmEmission(K_data::DataFrame)
    
    I_135nm =  K_data[!, @sprintf("r%i: O_5s -> O", ix_ref+253)]
    I_135nm += K_data[!, @sprintf("r%i: O_5s -> O", ix_ref+254)]
    I_135nm += K_data[!, @sprintf("r%i: e + O2 -> e + 2O", ix_ref+4)]
    I_135nm += K_data[!, @sprintf("r%i: e + O2_a1Ag -> e + 2O", ix_ref+5)]
    I_135nm += K_data[!, @sprintf("r%i: e + O2_b1Su -> e + 2O", ix_ref+6)]
    #I_135nm += K_data[!, @sprintf("r%i: e + O2_a1Ag_v -> e + 2O", ix_ref+8)]
    #I_135nm += K_data[!, @sprintf("r%i: e + O2_b1Su_v -> e + 2O", ix_ref+9)]
    return I_135nm
end

end