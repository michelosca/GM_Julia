module EvaluateExpressions

using SharedData: e, K_to_eV, amu, kb
using SharedData: Species, System, SpeciesID


function ReplaceTempSymbolS!(expr::Expr)

    ReplaceSymbol!(expr, :Te,     :(temp[sID.electron]))
    ReplaceSymbol!(expr, :TAr,    :(temp[sID.Ar]))
    ReplaceSymbol!(expr, :TArIon, :(temp[sID.Ar_Ion]))
    ReplaceSymbol!(expr, :TO2,    :(temp[sID.O2]))
    ReplaceSymbol!(expr, :TO,     :(temp[sID.O]))
end

function ReplaceDensSymbolS!(expr::Expr)

    ReplaceSymbol!(expr, :n_O2,     :(dens[sID.O2]))

end

function ReplaceSpeciesSymbolS!(expr::Expr)
    # Species parameters

    ReplaceSymbol!(expr, :uB_0,        :(species_list[sID.O].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O,       :(species_list[sID.O].v_thermal) )
    ReplaceSymbol!(expr, :h_R_O,       :(species_list[sID.O].h_R) )
    ReplaceSymbol!(expr, :h_L_O,       :(species_list[sID.O].h_L) )
    ReplaceSymbol!(expr, :D_O,         :(species_list[sID.O].D) )
    ReplaceSymbol!(expr, :gamma_O,    :(species_list[sID.O].gamma) )
    
    ReplaceSymbol!(expr, :uB_02,       :(species_list[sID.O].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2,      :(species_list[sID.O2].v_thermal) )
    ReplaceSymbol!(expr, :h_R_O2,      :(species_list[sID.O2].h_R) )
    ReplaceSymbol!(expr, :h_L_O2,      :(species_list[sID.O2].h_L) )
    ReplaceSymbol!(expr, :D_O2,        :(species_list[sID.O2].D) )
    ReplaceSymbol!(expr, :gamma_O2,    :(species_list[sID.O2].gamma) )
end


function ReplaceSystemSymbolS!(expr::Expr)
    # System parameters
    ReplaceSymbol!(expr, :R,         :(system.radius) )
    ReplaceSymbol!(expr, :L,         :(system.l) )
    ReplaceSymbol!(expr, :A,         :(system.A) )
    ReplaceSymbol!(expr, :V,         :(system.V) )
    ReplaceSymbol!(expr, :Lambda,    :(system.Lambda) )
end


function ReplaceConstantSymbolS!(expr::Expr)
    ReplaceSymbol!(expr, :m_Ar,      Expr(:call,:*,amu, 40 ))
    ReplaceSymbol!(expr, :m_O,       Expr(:call,:*,amu, 16 ))
    ReplaceSymbol!(expr, :m_O2,      Expr(:call,:*,amu, 32 ))
    ReplaceSymbol!(expr, :Te_eV,     Expr(:call,:*,:Te, K_to_eV))
    ReplaceSymbol!(expr, :TArIon_eV, Expr(:call,:*,:TArIon, K_to_eV))
    ReplaceSymbol!(expr, :kb,        kb)
    ReplaceSymbol!(expr, :e,         e)
    ReplaceSymbol!(expr, :pi,        pi)
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