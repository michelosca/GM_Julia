#push!(LOAD_PATH, pwd())
using SharedData
using PrintModule
using GenerateODEs
using ReadInput
using Plots

using SharedData: e, mp

global dens = Float64[]
global temp = Float64[]

function run_GM()
    # Reads data from the input.deck file
    errcode = ReadInput.SetupInputData()
    if (errcode != 0)
        return errcode
    end

    # Print species and reaction lists to terminal
    PrintModule.PrintSpeciesList()
    PrintModule.PrintReactionList()

    dens_funct_list = GenerateODEs.GenerateDensRateFunctionList()
    temp_funct_list = GenerateODEs.GenerateTempRateFunctionList()

    # Initial conditions
    # Density
    n_e = 1.e14
    n_Ar = 1.e20
    n_ArIon = n_e
    dens = [n_e, n_Ar, n_ArIon]

    # Temp
    T_e = 3 * SharedData.e / SharedData.kb
    T_Ar = 300
    T_ArIon = T_Ar
    temp = [T_e, T_Ar, T_ArIon]

    GenerateODEs.FirstUpdateGlobalSymbols(temp::Vector{Float64})
    #TestRateCoefficients()
    TestDensEquations(dens, temp, dens_funct_list)
    #TestTempEquations()
end

function TestDensEquations(dens, temp, dens_eq)
    function ne(dens, temp)
        Te = temp[SharedData.s_electron_id] 
        Te_eV = Te * SharedData.K_to_eV
        n_e = dens[SharedData.s_electron_id]
        n_Ar = dens[SharedData.s_Ar_id]
        n_ArIon = dens[SharedData.s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        m_Ar = 4 * SharedData.mp
        K_elastic = 1.0e-18 * sqrt(((8 * e * Te_eV) / pi) / m_Ar)
        n_gain = n_e * n_Ar * K_ionization
        n_loss = n_e * n_ArIon * K_recombination
        flux_e = n_e * pi * SharedData.kb * Te / 
            (GenerateODEs.m_Ar * n_Ar * K_elastic * SharedData.l) 
        n_flux = SharedData.A / SharedData.V * flux_e
        f = n_gain - n_loss - n_flux
        return f
    end

    function nArIon(dens, temp)
        Te = temp[SharedData.s_electron_id] 
        Te_eV = Te * SharedData.K_to_eV
        n_e = dens[SharedData.s_electron_id]
        n_Ar = dens[SharedData.s_Ar_id]
        n_ArIon = dens[SharedData.s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        m_Ar = 4 * SharedData.mp
        K_elastic = 1.0e-18 * sqrt(((8 * e * Te_eV) / pi) / m_Ar)
        n_gain = n_e * n_Ar * K_ionization
        n_loss = n_e * n_ArIon * K_recombination
        flux_ArIon = n_e * pi * SharedData.kb * Te / 
            (GenerateODEs.m_Ar * n_Ar * K_elastic * SharedData.l) 
        n_flux = SharedData.A / SharedData.V * flux_ArIon
        f = n_gain - n_loss - n_flux
        return f
    end

    function nAr(dens, temp)
        Te = temp[SharedData.s_electron_id] 
        Te_eV = Te * SharedData.K_to_eV
        n_e = dens[SharedData.s_electron_id]
        n_Ar = dens[SharedData.s_Ar_id]
        n_ArIon = dens[SharedData.s_ArIon_id]
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
    ne_sim_eq = dens_eq[SharedData.s_electron_id]
    n_th = zeros(n)
    f_vals = zeros(n)
    # ArIon equation
    ArIon_sim_eq = dens_eq[SharedData.s_ArIon_id]
    nArIon_th = zeros(n)
    fArIon_vals = zeros(n)
    # Ar equation
    Ar_sim_eq = dens_eq[SharedData.s_Ar_id]
    nAr_th = zeros(n)
    fAr_vals = zeros(n)
    for i in 1:n
        temp[SharedData.s_electron_id] = T_vec[i] / SharedData.K_to_eV
        GenerateODEs.UpdateGlobalSymbols(temp)
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

end

function TestTempEquations()
    return 0
end
    
function TestRateCoefficients()
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
        m_Ar = 4 * SharedData.mp
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
    global m_Ar = 4 * mp
    for i in 1:length(SharedData.reaction_list)
        K_val = SharedData.reaction_list[i].rate_coefficient
        K_sim_vals = zeros(n)
        for i in 1:n
            global Te_eV = T_vec[i]
            K_sim_vals[i] = eval(K_val)
        end
        plot!(p_list[i], T_vec, K_sim_vals, label="Sim k$i")
    end
    display(plot(p1,p2,p3,p4))
end