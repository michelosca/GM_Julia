#push!(LOAD_PATH, pwd())
using SharedData
using PrintModule
using GenerateODEs
using ReadInput
using Plots

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
    dens, temp = SetupInitialConditions()
    print("Initial density values", dens,"\n")
    print("Initial temperature values", temp,"\n")
    
    #TestRateCoefficients()
    #TestDensEquations(dens, temp, dens_funct_list)
    #TestTempEquations(dens, temp, temp_funct_list)
end

function SetupInitialConditions()
    dens = Float64[]
    temp = Float64[]
    for s in SharedData.species_list
        push!(dens, s.dens0)
        push!(temp, s.temp0)
    end
    GenerateODEs.FirstUpdateGlobalSymbols(temp)
    return dens, temp
end

# TEST FUNCTIONS
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
        K_elastic = 1.0e-18 * sqrt(((8 * SharedData.e * Te_eV) / pi) / m_Ar)
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
        K_elastic = 1.0e-18 * sqrt(((8 * SharedData.e * Te_eV) / pi) / m_Ar)
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


    n = 300
    dens_vec = LogRange(1.e13,1.e25,n)
    temp[SharedData.s_electron_id] = 3 / SharedData.K_to_eV
    GenerateODEs.UpdateGlobalSymbols(temp)
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
        dens[SharedData.s_Ar_id] = dens_vec[i]
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
    dens[SharedData.s_Ar_id] = 1.e20 

end

function TestTempEquations(dens, temp, temp_eq)
    function Te_eq_th(dens, temp)
        Te = temp[SharedData.s_electron_id]
        TAr = temp[SharedData.s_Ar_id]
        Te_eV = Te * SharedData.K_to_eV
        n_e = dens[SharedData.s_electron_id]
        n_Ar = dens[SharedData.s_Ar_id]
        n_ArIon = dens[SharedData.s_ArIon_id]
        K_ionization = 2.34e-14 * Te_eV ^ 0.59 * exp(-17.44 / Te_eV)
        K_recombination = 5.0e-39 * Te_eV ^ 4.5
        m_Ar = 4 * SharedData.mp
        K_elastic = 1.0e-18 * sqrt(((8 * SharedData.e * Te_eV) / pi) / m_Ar)
        K_elastic_e = 2.336e-14 * Te_eV^1.609 * exp(0.0618 * (log(Te_eV))^2 - 0.1171*(log(Te_eV))^3)

        T_elastic = -3*SharedData.me/m_Ar * K_elastic_e*n_e * n_Ar*SharedData.kb*(Te - TAr)
        T_gain = n_e * n_Ar * K_ionization
        T_loss = n_e * n_ArIon * K_recombination
        T_Ethreshold = -n_e * n_Ar * K_ionization * 15.76 * SharedData.e
        flux_ArIon = n_e * pi * SharedData.kb * Te / 
            (GenerateODEs.m_Ar * n_Ar * K_elastic * SharedData.l)
        V_sheath = 3/4 * (SharedData.drivI/SharedData.A)^2 /
            (SharedData.e * SharedData.eps0 * n_e * SharedData.drivOmega^2)
        T_flux = -SharedData.A / SharedData.V * flux_ArIon*(2.5*SharedData.kb*Te + SharedData.e * V_sheath)
        S_abs = SharedData.drivP / SharedData.V
        c = 1.5 * SharedData.kb * n_e
        c1 = -1.5*SharedData.kb*Te
        f = (S_abs + T_elastic + (T_gain - T_loss)*c1 + T_Ethreshold + T_flux) / c
        return f
    end

    LogRange(x1, x2, n) = [10^y for y in range(log10(x1), log10(x2), length=n)]
    n = 300
    T_vec = LogRange(0.1,100,n)

    # electron equation
    Te_sim_eq = temp_eq[SharedData.s_electron_id]
    Te_th = zeros(n)
    fTe_vals = zeros(n)
    for i in 1:n
        temp[SharedData.s_electron_id] = T_vec[i] / SharedData.K_to_eV
        GenerateODEs.UpdateGlobalSymbols(temp)
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
        k = 1.e-18 * sqrt(8*SharedData.e*Te_eV/pi/m_Ar)
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
    global m_Ar = 4 * SharedData.mp
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