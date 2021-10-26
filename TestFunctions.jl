module TestFunctions

using SharedData: K_to_eV, amu, e, kb, me, eps0
using GenerateODEs
using InputBlock_Species: s_electron_id, s_Ar_id, s_ArIon_id
using Plots

# Following test functions work with the following reaction reaction set
#  # reaction process       # threshold E    # rate coefficient                                                                  # description
#  e + Ar -> e + Ar;        0;               2.336e-14 * Te_eV^1.609 * exp(0.0618 * (log(Te_eV))^2 - 0.1171*(log(Te_eV))^3);           elastic
#  e + Ar -> e + e + Ar+;   15.76 eV;        2.34e-14 * Te_eV^0.59 * exp(-17.44/Te_eV);                                             ionization
#  e + Ar+ -> Ar;           0;               5e-39 * Te_eV^4.5;                                                                  recombination
#  Ar + Ar+ -> Ar + Ar+;    0;               1.e-18 * sqrt(8*e*Te_eV/pi/m_Ar);                                                      elastic


###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - TestDensEquations
# - TestTempEquations
# - TestRateCoefficients


function TestDensEquations(dens, temp, dens_eq, system)
    function ne(dens, temp)
        Te = temp[s_electron_id] 
        Te_eV = Te * K_to_eV
        n_e = dens[s_electron_id]
        n_Ar = dens[s_Ar_id]
        n_ArIon = dens[s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        m_Ar = 4 * amu
        K_elastic = 1.0e-18 * sqrt(((8 * e * Te_eV) / pi) / m_Ar)
        n_gain = n_e * n_Ar * K_ionization
        n_loss = n_e * n_ArIon * K_recombination
        flux_e = n_e * pi * kb * Te / 
            (m_Ar * n_Ar * K_elastic * system.l) 
        n_flux = system.A / system.V * flux_e
        f = n_gain - n_loss - n_flux
        return f
    end

    function nArIon(dens, temp)
        Te = temp[s_electron_id] 
        Te_eV = Te * K_to_eV
        n_e = dens[s_electron_id]
        n_Ar = dens[s_Ar_id]
        n_ArIon = dens[s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        m_Ar = 4 * amu
        K_elastic = 1.0e-18 * sqrt(((8 * e * Te_eV) / pi) / m_Ar)
        n_gain = n_e * n_Ar * K_ionization
        n_loss = n_e * n_ArIon * K_recombination
        flux_ArIon = n_e * pi * kb * Te / 
            (m_Ar * n_Ar * K_elastic * system.l) 
        n_flux = system.A / system.V * flux_ArIon
        f = n_gain - n_loss - n_flux
        return f
    end

    function nAr(dens, temp)
        Te = temp[s_electron_id] 
        Te_eV = Te * K_to_eV
        n_e = dens[s_electron_id]
        n_Ar = dens[s_Ar_id]
        n_ArIon = dens[s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        n_loss = n_e * n_Ar * K_ionization
        n_gain = n_e * n_ArIon * K_recombination
        f = n_gain - n_loss
        return f
    end

    LogRange(x1, x2, n) = [10^y for y in range(log10(x1), log10(x2), length=n)]
    n = 300
    T_vec = LogRange(0.1,100,n)

    #Electron equation
    ne_sim_eq = dens_eq[s_electron_id]
    n_th = zeros(n)
    f_vals = zeros(n)
    # ArIon equation
    ArIon_sim_eq = dens_eq[s_ArIon_id]
    nArIon_th = zeros(n)
    fArIon_vals = zeros(n)
    # Ar equation
    Ar_sim_eq = dens_eq[s_Ar_id]
    nAr_th = zeros(n)
    fAr_vals = zeros(n)
    for i in 1:n
        temp[s_electron_id] = T_vec[i] / K_to_eV
        #Electrons###
        # Theory
        n_th[i] = abs(ne(dens, temp))
        # Sim
        f_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, ne_sim_eq))

        # ArIons
        # Theory
        nArIon_th[i] = abs(nArIon(dens, temp))
        # Sim
        fArIon_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, ArIon_sim_eq))

        # Ar
        # Theory
        nAr_th[i] = abs(nAr(dens, temp))
        # Sim
        fAr_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, Ar_sim_eq))
    end
    p1 = plot(T_vec, n_th, xaxis=:log, yaxis=:log, label="Electron dens Th")
    plot!(p1, T_vec, f_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="Electron dens Sim")
    display(plot(p1))
    p2 = plot(T_vec, nArIon_th, xaxis=:log, yaxis=:log, label="ArIon dens Th")
    plot!(p2, T_vec, fArIon_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="ArIon dens Sim")
    display(plot(p2))
    p3 = plot(T_vec, nAr_th, xaxis=:log, yaxis=:log, label="Ar dens Th")
    plot!(p3, T_vec, fAr_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="Ar dens Sim")
    display(plot(p3))


    n = 300
    dens_vec = LogRange(1.e13,1.e25,n)
    temp[s_electron_id] = 3 / K_to_eV
    #Electron equation
    ne_sim_eq = dens_eq[s_electron_id]
    n_th = zeros(n)
    f_vals = zeros(n)
    # ArIon equation
    ArIon_sim_eq = dens_eq[s_ArIon_id]
    nArIon_th = zeros(n)
    fArIon_vals = zeros(n)
    # Ar equation
    Ar_sim_eq = dens_eq[s_Ar_id]
    nAr_th = zeros(n)
    fAr_vals = zeros(n)
    for i in 1:n
        dens[s_Ar_id] = dens_vec[i]
        #Electrons###
        # Theory
        n_th[i] = abs(ne(dens, temp))
        # Sim
        f_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, ne_sim_eq))

        # ArIons
        # Theory
        nArIon_th[i] = abs(nArIon(dens, temp))
        # Sim
        fArIon_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, ArIon_sim_eq))

        # Ar
        # Theory
        nAr_th[i] = abs(nAr(dens, temp))
        # Sim
        fAr_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, Ar_sim_eq))
    end
    p1 = plot(dens_vec, n_th, xaxis=:log, yaxis=:log, label="Electron dens Th")
    plot!(p1, dens_vec, f_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="Electron dens Sim")
    display(plot(p1))
    p2 = plot(dens_vec, nArIon_th, xaxis=:log, yaxis=:log, label="ArIon dens Th")
    plot!(p2, dens_vec, fArIon_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="ArIon dens Sim")
    display(plot(p2))
    p3 = plot(dens_vec, nAr_th, xaxis=:log, yaxis=:log, label="Ar dens Th")
    plot!(p3, dens_vec, fAr_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="Ar dens Sim")
    display(plot(p3))
    dens[s_Ar_id] = 1.e20 

end

function TestTempEquations(dens, temp, temp_eq, system)
    function Te_eq_th(dens, temp)
        Te = temp[s_electron_id]
        TAr = temp[s_Ar_id]
        Te_eV = Te * K_to_eV
        n_e = dens[s_electron_id]
        n_Ar = dens[s_Ar_id]
        n_ArIon = dens[s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        m_Ar = 4 * amu
        K_elastic = 1.0e-18 * sqrt(((8 * e * Te_eV) / pi) / m_Ar)
        K_elastic_e = 2.336e-14 * Te_eV^1.609 * exp(0.0618 * (log(Te_eV))^2 - 0.1171*(log(Te_eV))^3)

        T_elastic = -3*me/m_Ar * K_elastic_e*n_e * n_Ar*kb*(Te - TAr)
        T_gain = n_e * n_Ar * K_ionization
        T_loss = n_e * n_ArIon * K_recombination
        T_Ethreshold = -n_e * n_Ar * K_ionization * 15.76 * e
        flux_ArIon = n_e * pi * kb * Te / 
            (m_Ar * n_Ar * K_elastic * system.l)
        V_sheath = 3/4 * (system.drivI/system.A)^2 /
            (e * eps0 * n_e * system.drivOmega^2)
        T_flux = -system.A / system.V * flux_ArIon*(2.5*kb*Te + e * V_sheath)
        S_abs = system.drivP / system.V
        c = 1.5 * kb * n_e
        c1 = -1.5*kb*Te
        f = (S_abs + T_elastic + (T_gain - T_loss)*c1 + T_Ethreshold + T_flux) / c
        return f
    end

    LogRange(x1, x2, n) = [10^y for y in range(log10(x1), log10(x2), length=n)]
    n = 300
    T_vec = LogRange(0.1,100,n)

    # electron equation
    Te_sim_eq = temp_eq[s_electron_id]
    Te_th = zeros(n)
    fTe_vals = zeros(n)
    for i in 1:n
        temp[s_electron_id] = T_vec[i] / K_to_eV
        #Electrons###
        # Theory
        Te_th[i] = abs(Te_eq_th(dens, temp))
        # Sim
        fTe_vals[i] = abs(GenerateODEs.GatherListOfFunctions(dens, temp, Te_sim_eq))
    end
    p1 = plot(T_vec, Te_th, xaxis=:log, yaxis=:log, label="Electron temp Th")
    plot!(p1, T_vec, fTe_vals, xaxis=:log, yaxis=:log, line = (:dot, 4),label="Electron temp Sim")
    display(plot(p1))
end
    
function TestRateCoefficients(temp, reaction_list)
    # Tests function
    function K1(Te_eV)
        k = 2.336e-14 * Te_eV^1.609 * exp(0.0618 * (log(Te_eV))^2 - 0.1171*(log(Te_eV))^3)
        return k
    end
    function K2(Te_eV)
        k =  2.34e-14 * Te_eV^0.59 * exp(-17.44/Te_eV)
        return k
    end
    function K3(Te_eV)
        k = 5e-39 * Te_eV^4.5
        return k
    end
    function K4(Te_eV)
        m_Ar = 4 * amu
        k = 1.e-18 * sqrt(8*e*Te_eV/pi/m_Ar)
        return k
    end
    LogRange(x1, x2, n) = [10^y for y in range(log10(x1), log10(x2), length=n)]
    n = 300
    T_vec = LogRange(0.1,10000,n)
    p1 = plot(T_vec, map(K1, T_vec), xaxis=:log, yaxis=:log, label="Input k1")
    p2 = plot(T_vec, map(K2, T_vec), xaxis=:log, yaxis=:log, label="Input k2")
    p3 = plot(T_vec, map(K3, T_vec), xaxis=:log, yaxis=:log, label="Input k3")
    p4 = plot(T_vec, map(K4, T_vec), xaxis=:log, yaxis=:log, label="Input k4")
    p_list = [p1,p2,p3,p4]

    # Sim function
    for i in 1:length(reaction_list)
        K_sim_vals = zeros(n)
        for j in 1:n
            temp[1] = T_vec[j] * 1.16e4
            K_sim_vals[j] = reaction_list[i].rate_coefficient(temp)
        end
        plot!(p_list[i], T_vec, K_sim_vals, line = (:dot, 4), label="Sim k$i")
    end
    display(plot(p1,p2,p3,p4))
end

end