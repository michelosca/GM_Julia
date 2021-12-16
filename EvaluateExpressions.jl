module EvaluateExpressions

using SharedData: e, K_to_eV, amu, kb
using SharedData: Species, System, SpeciesID

function ReplaceExpressionValues(expr::Union{Float64, Expr}, temp::Vector{Float64},
    species_list::Vector{Species}, system::System, sID::SpeciesID)

    copy_expr = copy(expr)
    if typeof(copy_expr)==Expr
        ReplaceSymbolValues!(copy_expr, temp, species_list, system, sID)
    end

    value = eval(copy_expr)
    return value

end


function ReplaceTempSymbols!(expr::Expr)

    ReplaceSymbol!(expr, :T_e,      :(temp[sID.electron]))
    ReplaceSymbol!(expr, :T_Ar,     :(temp[sID.Ar]))
    ReplaceSymbol!(expr, :T_Ar_Ion, :(temp[sID.Ar_Ion]))
    ReplaceSymbol!(expr, :T_O,      :(temp[sID.O]))
    ReplaceSymbol!(expr, :T_O2,     :(temp[sID.O2]))
    ReplaceSymbol!(expr, :T_O3,     :(temp[sID.O3]))
end


function ReplaceSymbolValues!(expr::Expr, temp::Vector{Float64},
    species_list::Vector{Species}, system::System,
    sID::SpeciesID)

    # electrons
    ReplaceSymbol!(expr, :T_e,      temp[sID.electron])

    # Ar neutrals
    if sID.Ar != 0
        ReplaceSymbol!(expr, :T_Ar,     temp[sID.Ar])
        ReplaceSymbol!(expr, :uB_Ar,    species_list[sID.Ar].v_Bohm)
        ReplaceSymbol!(expr, :vth_Ar,   species_list[sID.Ar].v_thermal)
        ReplaceSymbol!(expr, :h_R_Ar,   species_list[sID.Ar].h_R)
        ReplaceSymbol!(expr, :h_L_Ar,   species_list[sID.Ar].h_L)
        ReplaceSymbol!(expr, :D_Ar,     species_list[sID.Ar].D)
        ReplaceSymbol!(expr, :gamma_Ar, species_list[sID.Ar].gamma)
    end

    # Ar ions
    if sID.Ar_Ion != 0
        ReplaceSymbol!(expr, :T_Ar_Ion, temp[sID.Ar_Ion])
    end

    # O2 neutrals
    if sID.O2 != 0
        ReplaceSymbol!(expr, :T_O2,     temp[sID.O2])
        ReplaceSymbol!(expr, :uB_O2,    species_list[sID.O2].v_Bohm)
        ReplaceSymbol!(expr, :vth_O2,   species_list[sID.O2].v_thermal)
        ReplaceSymbol!(expr, :h_R_O2,   species_list[sID.O2].h_R)
        ReplaceSymbol!(expr, :h_L_O2,   species_list[sID.O2].h_L)
        ReplaceSymbol!(expr, :D_O2,     species_list[sID.O2].D)
        ReplaceSymbol!(expr, :gamma_O2, species_list[sID.O2].gamma)
    end

    # O neutrals
    if sID.O != 0
        ReplaceSymbol!(expr, :T_O,      temp[sID.O])
        ReplaceSymbol!(expr, :uB_O,     species_list[sID.O].v_Bohm)
        ReplaceSymbol!(expr, :vth_O,    species_list[sID.O].v_thermal)
        ReplaceSymbol!(expr, :h_R_O,    species_list[sID.O].h_R)
        ReplaceSymbol!(expr, :h_L_O,    species_list[sID.O].h_L)
        ReplaceSymbol!(expr, :D_O,      species_list[sID.O].D)
        ReplaceSymbol!(expr, :gamma_O,  species_list[sID.O].gamma)
    end

    # O3 neutrals
    if sID.O3 != 0
        ReplaceSymbol!(expr, :T_O3,     temp[sID.O3])
        ReplaceSymbol!(expr, :uB_O3,    species_list[sID.O3].v_Bohm)
        ReplaceSymbol!(expr, :vth_O3,   species_list[sID.O3].v_thermal)
        ReplaceSymbol!(expr, :h_R_O3,   species_list[sID.O3].h_R)
        ReplaceSymbol!(expr, :h_L_O3,   species_list[sID.O3].h_L)
        ReplaceSymbol!(expr, :D_O3,     species_list[sID.O3].D)
        ReplaceSymbol!(expr, :gamma_O3, species_list[sID.O3].gamma)
    end

    # System values
    ReplaceSymbol!(expr, :R,         system.radius)
    ReplaceSymbol!(expr, :L,         system.l)
    ReplaceSymbol!(expr, :A,         system.A)
    ReplaceSymbol!(expr, :V,         system.V)
    ReplaceSymbol!(expr, :Lambda,    system.Lambda)
end


function ReplaceDensSymbols!(expr::Expr)
    ReplaceSymbol!(expr, :n_O2,     :(dens[sID.O2]))
end


function ReplaceSpeciesSymbols!(expr::Expr)

    # Species parameters
    ReplaceSymbol!(expr, :uB_O,        :(species_list[sID.O].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O,       :(species_list[sID.O].v_thermal) )
    ReplaceSymbol!(expr, :h_R_O,       :(species_list[sID.O].h_R) )
    ReplaceSymbol!(expr, :h_L_O,       :(species_list[sID.O].h_L) )
    ReplaceSymbol!(expr, :D_O,         :(species_list[sID.O].D) )
    ReplaceSymbol!(expr, :gamma_O,     :(species_list[sID.O].gamma) )
    
    ReplaceSymbol!(expr, :uB_O2,       :(species_list[sID.O2].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2,      :(species_list[sID.O2].v_thermal) )
    ReplaceSymbol!(expr, :h_R_O2,      :(species_list[sID.O2].h_R) )
    ReplaceSymbol!(expr, :h_L_O2,      :(species_list[sID.O2].h_L) )
    ReplaceSymbol!(expr, :D_O2,        :(species_list[sID.O2].D) )
    ReplaceSymbol!(expr, :gamma_O2,    :(species_list[sID.O2].gamma) )

    ReplaceSymbol!(expr, :uB_O3,       :(species_list[sID.O3].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O3,      :(species_list[sID.O3].v_thermal) )
    ReplaceSymbol!(expr, :h_R_O3,      :(species_list[sID.O3].h_R) )
    ReplaceSymbol!(expr, :h_L_O3,      :(species_list[sID.O3].h_L) )
    ReplaceSymbol!(expr, :D_O3,        :(species_list[sID.O3].D) )
    ReplaceSymbol!(expr, :gamma_O3,    :(species_list[sID.O3].gamma) )

    ReplaceSymbol!(expr, :uB_Ar,       :(species_list[sID.Ar].v_Bohm) )
    ReplaceSymbol!(expr, :vth_Ar,      :(species_list[sID.Ar].v_thermal) )
    ReplaceSymbol!(expr, :h_R_Ar,      :(species_list[sID.Ar].h_R) )
    ReplaceSymbol!(expr, :h_L_Ar,      :(species_list[sID.Ar].h_L) )
    ReplaceSymbol!(expr, :D_Ar,        :(species_list[sID.Ar].D) )
    ReplaceSymbol!(expr, :gamma_Ar,    :(species_list[sID.Ar].gamma) )
end


function ReplaceSystemSymbols!(expr::Expr)
    # System parameters
    ReplaceSymbol!(expr, :R,         :(system.radius) )
    ReplaceSymbol!(expr, :L,         :(system.l) )
    ReplaceSymbol!(expr, :A,         :(system.A) )
    ReplaceSymbol!(expr, :V,         :(system.V) )
    ReplaceSymbol!(expr, :Lambda,    :(system.Lambda) )
end


function ReplaceConstantValues!(expr::Expr)
    ReplaceSymbol!(expr, :m_Ar,        Expr(:call,:*,amu, 40 ))
    ReplaceSymbol!(expr, :m_O,         Expr(:call,:*,amu, 16 ))
    ReplaceSymbol!(expr, :m_O2,        Expr(:call,:*,amu, 32 ))
    ReplaceSymbol!(expr, :m_O3,        Expr(:call,:*,amu, 48 ))
    ReplaceSymbol!(expr, :m_O4,        Expr(:call,:*,amu, 64 ))
    ReplaceSymbol!(expr, :T_e_eV,      Expr(:call,:*,:T_e, K_to_eV))
    ReplaceSymbol!(expr, :T_Ar_Ion_eV, Expr(:call,:*,:T_Ar_Ion, K_to_eV))
    ReplaceSymbol!(expr, :kb,          kb)
    ReplaceSymbol!(expr, :e,           e)
    ReplaceSymbol!(expr, :pi,          pi)
end


function ReplaceSymbol!(expression::Expr, old, new)
    n = length(expression.args)
    for i in 1:n
        arg = expression.args[i]
        if (typeof(arg) == Expr)
            ReplaceSymbol!(arg, old, new)
        elseif (arg == old)
            expression.args[i] = new
        end
    end
end


end