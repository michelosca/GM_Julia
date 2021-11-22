module EvaluateExpressions

using SharedData: e, K_to_eV, amu, kb
using SharedData: Species, System, SpeciesID


function EvaluateExpression(expr, sID::SpeciesID, temp::Vector{Float64},
    species::Species)

    new_expr = copy(expr)
    if !(typeof(new_expr) == Float64)
        ReplaceTempSymbolS!(new_expr, sID, temp)
        ReplaceSpeciesSymbolS!(new_expr, species)
    end

    return new_expr 
end


function ReplaceTempSymbolS!(expr::Expr, sID::SpeciesID, temp::Vector{Float64})

    ReplaceSymbol!(expr, :Te,     temp[sID.electron])
    if sID.Ar != 0
        ReplaceSymbol!(expr, :TAr,     temp[sID.Ar])
    end
    if sID.Ar_Ion != 0
        ReplaceSymbol!(expr, :TArIon, temp[sID.Ar_Ion])
    end
    if sID.O2 != 0
        ReplaceSymbol!(expr, :TO2,    temp[sID.O2])
    end
    if sID.O != 0
        ReplaceSymbol!(expr, :TO,     temp[sID.O])
    end
end


function ReplaceSpeciesSymbolS!(expr::Expr, species::Species)
    # Species parameters
    ReplaceSymbol!(expr, :h_R,       species.h_R )
    ReplaceSymbol!(expr, :h_L,       species.h_L )
    ReplaceSymbol!(expr, :uB,        species.v_Bohm )
    ReplaceSymbol!(expr, :vth,       species.v_thermal )
    ReplaceSymbol!(expr, :Lambda,    species.Lambda )
    ReplaceSymbol!(expr, :D,         species.D )
    ReplaceSymbol!(expr, :gamma,     species.gamma )
end


function ReplaceSystemSymbolS!(expr::Expr, system::System)
    # System parameters
    ReplaceSymbol!(expr, :R,         system.radius )
    ReplaceSymbol!(expr, :L,         system.l )
    ReplaceSymbol!(expr, :A,         system.A )
    ReplaceSymbol!(expr, :V,         system.V )
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