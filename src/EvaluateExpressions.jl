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
    ReplaceSymbol!(expr, :T_O_1d,   :(temp[sID.O_1d]))
    ReplaceSymbol!(expr, :T_O_1s,   :(temp[sID.O_1s]))
    ReplaceSymbol!(expr, :T_O_3s,   :(temp[sID.O_3s]))
    ReplaceSymbol!(expr, :T_O_5s,   :(temp[sID.O_5s]))
    ReplaceSymbol!(expr, :T_O_3p,   :(temp[sID.O_3s]))
    ReplaceSymbol!(expr, :T_O_5p,   :(temp[sID.O_5s]))
    ReplaceSymbol!(expr, :T_O2,     :(temp[sID.O2]))
    ReplaceSymbol!(expr, :T_O3,     :(temp[sID.O3]))
    ReplaceSymbol!(expr, :T_O4,     :(temp[sID.O4]))
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
        ReplaceSymbol!(expr, :D_Ar,     species_list[sID.Ar].D)
        ReplaceSymbol!(expr, :gamma_Ar, species_list[sID.Ar].gamma)
    end
    if sID.Ar_m != 0
        ReplaceSymbol!(expr, :T_Ar_m,     temp[sID.Ar_m])
        ReplaceSymbol!(expr, :uB_Ar_m,    species_list[sID.Ar_m].v_Bohm)
        ReplaceSymbol!(expr, :vth_Ar_m,   species_list[sID.Ar_m].v_thermal)
        ReplaceSymbol!(expr, :D_Ar_m,     species_list[sID.Ar_m].D)
        ReplaceSymbol!(expr, :gamma_Ar_m, species_list[sID.Ar_m].gamma)
    end
    if sID.Ar_r != 0
        ReplaceSymbol!(expr, :T_Ar_r,     temp[sID.Ar_r])
        ReplaceSymbol!(expr, :uB_Ar_r,    species_list[sID.Ar_r].v_Bohm)
        ReplaceSymbol!(expr, :vth_Ar_r,   species_list[sID.Ar_r].v_thermal)
        ReplaceSymbol!(expr, :D_Ar_r,     species_list[sID.Ar_r].D)
        ReplaceSymbol!(expr, :gamma_Ar_r, species_list[sID.Ar_r].gamma)
    end
    if sID.Ar_4p != 0
        ReplaceSymbol!(expr, :T_Ar_4p,     temp[sID.Ar_4p])
        ReplaceSymbol!(expr, :uB_Ar_4p,    species_list[sID.Ar_4p].v_Bohm)
        ReplaceSymbol!(expr, :vth_Ar_4p,   species_list[sID.Ar_4p].v_thermal)
        ReplaceSymbol!(expr, :D_Ar_4p,     species_list[sID.Ar_4p].D)
        ReplaceSymbol!(expr, :gamma_Ar_4p, species_list[sID.Ar_4p].gamma)
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
        ReplaceSymbol!(expr, :D_O2,     species_list[sID.O2].D)
        ReplaceSymbol!(expr, :gamma_O2, species_list[sID.O2].gamma)
    end
    if sID.O2_v != 0
        ReplaceSymbol!(expr, :T_O2_v,     temp[sID.O2_v])
        ReplaceSymbol!(expr, :uB_O2_v,    species_list[sID.O2_v].v_Bohm)
        ReplaceSymbol!(expr, :vth_O2_v,   species_list[sID.O2_v].v_thermal)
        ReplaceSymbol!(expr, :D_O2_v,     species_list[sID.O2_v].D)
        ReplaceSymbol!(expr, :gamma_O2_v, species_list[sID.O2_v].gamma)
    end
    if sID.O2_a1Ag != 0
        ReplaceSymbol!(expr, :T_O2_a1Ag,     temp[sID.O2_a1Ag])
        ReplaceSymbol!(expr, :uB_O2_a1Ag,    species_list[sID.O2_a1Ag].v_Bohm)
        ReplaceSymbol!(expr, :vth_O2_a1Ag,   species_list[sID.O2_a1Ag].v_thermal)
        ReplaceSymbol!(expr, :D_O2_a1Ag,     species_list[sID.O2_a1Ag].D)
        ReplaceSymbol!(expr, :gamma_O2_a1Ag, species_list[sID.O2_a1Ag].gamma)
    end
    if sID.O2_b1Su != 0
        ReplaceSymbol!(expr, :T_O2_b1Su,     temp[sID.O2_b1Su])
        ReplaceSymbol!(expr, :uB_O2_b1Su,    species_list[sID.O2_b1Su].v_Bohm)
        ReplaceSymbol!(expr, :vth_O2_b1Su,   species_list[sID.O2_b1Su].v_thermal)
        ReplaceSymbol!(expr, :D_O2_b1Su,     species_list[sID.O2_b1Su].D)
        ReplaceSymbol!(expr, :gamma_O2_b1Su, species_list[sID.O2_b1Su].gamma)
    end
    if sID.O2_a1Ag_v != 0
        ReplaceSymbol!(expr, :T_O2_a1Ag_v,     temp[sID.O2_a1Ag_v])
        ReplaceSymbol!(expr, :uB_O2_a1Ag_v,    species_list[sID.O2_a1Ag_v].v_Bohm)
        ReplaceSymbol!(expr, :vth_O2_a1Ag_v,   species_list[sID.O2_a1Ag_v].v_thermal)
        ReplaceSymbol!(expr, :D_O2_a1Ag_v,     species_list[sID.O2_a1Ag_v].D)
        ReplaceSymbol!(expr, :gamma_O2_a1Ag_v, species_list[sID.O2_a1Ag_v].gamma)
    end
    if sID.O2_b1Su_v != 0
        ReplaceSymbol!(expr, :T_O2_b1Su_v,     temp[sID.O2_b1Su_v])
        ReplaceSymbol!(expr, :uB_O2_b1Su_v,    species_list[sID.O2_b1Su_v].v_Bohm)
        ReplaceSymbol!(expr, :vth_O2_b1Su_v,   species_list[sID.O2_b1Su_v].v_thermal)
        ReplaceSymbol!(expr, :D_O2_b1Su_v,     species_list[sID.O2_b1Su_v].D)
        ReplaceSymbol!(expr, :gamma_O2_b1Su_v, species_list[sID.O2_b1Su_v].gamma)
    end

    # O neutrals
    if sID.O != 0
        ReplaceSymbol!(expr, :T_O,      temp[sID.O])
        ReplaceSymbol!(expr, :uB_O,     species_list[sID.O].v_Bohm)
        ReplaceSymbol!(expr, :vth_O,    species_list[sID.O].v_thermal)
        ReplaceSymbol!(expr, :D_O,      species_list[sID.O].D)
        ReplaceSymbol!(expr, :gamma_O,  species_list[sID.O].gamma)
    end
    if sID.O_1d != 0
        ReplaceSymbol!(expr, :T_O_1d,      temp[sID.O_1d])
        ReplaceSymbol!(expr, :uB_O_1d,     species_list[sID.O_1d].v_Bohm)
        ReplaceSymbol!(expr, :vth_O_1d,    species_list[sID.O_1d].v_thermal)
        ReplaceSymbol!(expr, :D_O_1d,      species_list[sID.O_1d].D)
        ReplaceSymbol!(expr, :gamma_O_1d,  species_list[sID.O_1d].gamma)
    end
    if sID.O_1s != 0
        ReplaceSymbol!(expr, :T_O_1s,      temp[sID.O_1s])
        ReplaceSymbol!(expr, :uB_O_1s,     species_list[sID.O_1s].v_Bohm)
        ReplaceSymbol!(expr, :vth_O_1s,    species_list[sID.O_1s].v_thermal)
        ReplaceSymbol!(expr, :D_O_1s,      species_list[sID.O_1s].D)
        ReplaceSymbol!(expr, :gamma_O_1s,  species_list[sID.O_1s].gamma)
    end
    if sID.O_3s != 0
        ReplaceSymbol!(expr, :T_O_3s,      temp[sID.O_3s])
        ReplaceSymbol!(expr, :uB_O_3s,     species_list[sID.O_3s].v_Bohm)
        ReplaceSymbol!(expr, :vth_O_3s,    species_list[sID.O_3s].v_thermal)
        ReplaceSymbol!(expr, :D_O_3s,      species_list[sID.O_3s].D)
        ReplaceSymbol!(expr, :gamma_O_3s,  species_list[sID.O_3s].gamma)
    end
    if sID.O_5s != 0
        ReplaceSymbol!(expr, :T_O_5s,      temp[sID.O_5s])
        ReplaceSymbol!(expr, :uB_O_5s,     species_list[sID.O_5s].v_Bohm)
        ReplaceSymbol!(expr, :vth_O_5s,    species_list[sID.O_5s].v_thermal)
        ReplaceSymbol!(expr, :D_O_5s,      species_list[sID.O_5s].D)
        ReplaceSymbol!(expr, :gamma_O_5s,  species_list[sID.O_5s].gamma)
    end
    if sID.O_3p != 0
        ReplaceSymbol!(expr, :T_O_3p,      temp[sID.O_3p])
        ReplaceSymbol!(expr, :uB_O_3p,     species_list[sID.O_3p].v_Bohm)
        ReplaceSymbol!(expr, :vth_O_3p,    species_list[sID.O_3p].v_thermal)
        ReplaceSymbol!(expr, :D_O_3p,      species_list[sID.O_3p].D)
        ReplaceSymbol!(expr, :gamma_O_3p,  species_list[sID.O_3p].gamma)
    end
    if sID.O_5p != 0
        ReplaceSymbol!(expr, :T_O_5p,      temp[sID.O_5p])
        ReplaceSymbol!(expr, :uB_O_5p,     species_list[sID.O_5p].v_Bohm)
        ReplaceSymbol!(expr, :vth_O_5p,    species_list[sID.O_5p].v_thermal)
        ReplaceSymbol!(expr, :D_O_5p,      species_list[sID.O_5p].D)
        ReplaceSymbol!(expr, :gamma_O_5p,  species_list[sID.O_5p].gamma)
    end

    # O3 neutrals
    if sID.O3 != 0
        ReplaceSymbol!(expr, :T_O3,     temp[sID.O3])
        ReplaceSymbol!(expr, :uB_O3,    species_list[sID.O3].v_Bohm)
        ReplaceSymbol!(expr, :vth_O3,   species_list[sID.O3].v_thermal)
        ReplaceSymbol!(expr, :D_O3,     species_list[sID.O3].D)
        ReplaceSymbol!(expr, :gamma_O3, species_list[sID.O3].gamma)
    end
    if sID.O3_v != 0
        ReplaceSymbol!(expr, :T_O3_v,     temp[sID.O3_v])
        ReplaceSymbol!(expr, :uB_O3_v,    species_list[sID.O3_v].v_Bohm)
        ReplaceSymbol!(expr, :vth_O3_v,   species_list[sID.O3_v].v_thermal)
        ReplaceSymbol!(expr, :D_O3_v,     species_list[sID.O3_v].D)
        ReplaceSymbol!(expr, :gamma_O3_v, species_list[sID.O3_v].gamma)
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
    ReplaceSymbol!(expr, :n_O,     :(dens[sID.O]))
end


function ReplaceSpeciesSymbols!(expr::Expr)

    # Species parameters
    ReplaceSymbol!(expr, :uB_O,        :(species_list[sID.O].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O,       :(species_list[sID.O].v_thermal) )
    ReplaceSymbol!(expr, :D_O,         :(species_list[sID.O].D) )
    ReplaceSymbol!(expr, :gamma_O,     :(species_list[sID.O].gamma) )
    
    ReplaceSymbol!(expr, :uB_O_1d,        :(species_list[sID.O_1d].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O_1d,       :(species_list[sID.O_1d].v_thermal) )
    ReplaceSymbol!(expr, :D_O_1d,         :(species_list[sID.O_1d].D) )
    ReplaceSymbol!(expr, :gamma_O_1d,     :(species_list[sID.O_1d].gamma) )
    
    ReplaceSymbol!(expr, :uB_O_1s,        :(species_list[sID.O_1s].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O_1s,       :(species_list[sID.O_1s].v_thermal) )
    ReplaceSymbol!(expr, :D_O_1s,         :(species_list[sID.O_1s].D) )
    ReplaceSymbol!(expr, :gamma_O_1s,     :(species_list[sID.O_1s].gamma) )
    
    ReplaceSymbol!(expr, :uB_O_3s,        :(species_list[sID.O_3s].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O_3s,       :(species_list[sID.O_3s].v_thermal) )
    ReplaceSymbol!(expr, :D_O_3s,         :(species_list[sID.O_3s].D) )
    ReplaceSymbol!(expr, :gamma_O_3s,     :(species_list[sID.O_3s].gamma) )
    
    ReplaceSymbol!(expr, :uB_O_5s,        :(species_list[sID.O_5s].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O_5s,       :(species_list[sID.O_5s].v_thermal) )
    ReplaceSymbol!(expr, :D_O_5s,         :(species_list[sID.O_5s].D) )
    ReplaceSymbol!(expr, :gamma_O_5s,     :(species_list[sID.O_5s].gamma) )
    
    ReplaceSymbol!(expr, :uB_O_3p,        :(species_list[sID.O_3p].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O_3p,       :(species_list[sID.O_3p].v_thermal) )
    ReplaceSymbol!(expr, :D_O_3p,         :(species_list[sID.O_3p].D) )
    ReplaceSymbol!(expr, :gamma_O_3p,     :(species_list[sID.O_3p].gamma) )
    
    ReplaceSymbol!(expr, :uB_O_5p,        :(species_list[sID.O_5p].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O_5p,       :(species_list[sID.O_5p].v_thermal) )
    ReplaceSymbol!(expr, :D_O_5p,         :(species_list[sID.O_5p].D) )
    ReplaceSymbol!(expr, :gamma_O_5p,     :(species_list[sID.O_5p].gamma) )
    
    ReplaceSymbol!(expr, :uB_O2,       :(species_list[sID.O2].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2,      :(species_list[sID.O2].v_thermal) )
    ReplaceSymbol!(expr, :D_O2,        :(species_list[sID.O2].D) )
    ReplaceSymbol!(expr, :gamma_O2,    :(species_list[sID.O2].gamma) )

    ReplaceSymbol!(expr, :uB_O2_v,       :(species_list[sID.O2_v].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2_v,      :(species_list[sID.O2_v].v_thermal) )
    ReplaceSymbol!(expr, :D_O2_v,        :(species_list[sID.O2_v].D) )
    ReplaceSymbol!(expr, :gamma_O2_v,    :(species_list[sID.O2_v].gamma) )

    ReplaceSymbol!(expr, :uB_O2_a1Ag,       :(species_list[sID.O2_a1Ag].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2_a1Ag,      :(species_list[sID.O2_a1Ag].v_thermal) )
    ReplaceSymbol!(expr, :D_O2_a1Ag,        :(species_list[sID.O2_a1Ag].D) )
    ReplaceSymbol!(expr, :gamma_O2_a1Ag,    :(species_list[sID.O2_a1Ag].gamma) )

    ReplaceSymbol!(expr, :uB_O2_b1Su,       :(species_list[sID.O2_b1Su].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2_b1Su,      :(species_list[sID.O2_b1Su].v_thermal) )
    ReplaceSymbol!(expr, :D_O2_b1Su,        :(species_list[sID.O2_b1Su].D) )
    ReplaceSymbol!(expr, :gamma_O2_b1Su,    :(species_list[sID.O2_b1Su].gamma) )

    ReplaceSymbol!(expr, :uB_O2_a1Ag_v,       :(species_list[sID.O2_a1Ag_v].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2_a1Ag_v,      :(species_list[sID.O2_a1Ag_v].v_thermal) )
    ReplaceSymbol!(expr, :D_O2_a1Ag_v,        :(species_list[sID.O2_a1Ag_v].D) )
    ReplaceSymbol!(expr, :gamma_O2_a1Ag_v,    :(species_list[sID.O2_a1Ag_v].gamma) )

    ReplaceSymbol!(expr, :uB_O2_b1Su_v,       :(species_list[sID.O2_b1Su_v].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O2_b1Su_v,      :(species_list[sID.O2_b1Su_v].v_thermal) )
    ReplaceSymbol!(expr, :D_O2_b1Su_v,        :(species_list[sID.O2_b1Su_v].D) )
    ReplaceSymbol!(expr, :gamma_O2_b1Su_v,    :(species_list[sID.O2_b1Su_v].gamma) )

    ReplaceSymbol!(expr, :uB_O3,       :(species_list[sID.O3].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O3,      :(species_list[sID.O3].v_thermal) )
    ReplaceSymbol!(expr, :D_O3,        :(species_list[sID.O3].D) )
    ReplaceSymbol!(expr, :gamma_O3,    :(species_list[sID.O3].gamma) )

    ReplaceSymbol!(expr, :uB_O3_v,       :(species_list[sID.O3_v].v_Bohm) )
    ReplaceSymbol!(expr, :vth_O3_v,      :(species_list[sID.O3_v].v_thermal) )
    ReplaceSymbol!(expr, :D_O3_v,        :(species_list[sID.O3_v].D) )
    ReplaceSymbol!(expr, :gamma_O3_v,    :(species_list[sID.O3_v].gamma) )

    ReplaceSymbol!(expr, :uB_Ar,       :(species_list[sID.Ar].v_Bohm) )
    ReplaceSymbol!(expr, :vth_Ar,      :(species_list[sID.Ar].v_thermal) )
    ReplaceSymbol!(expr, :D_Ar,        :(species_list[sID.Ar].D) )
    ReplaceSymbol!(expr, :gamma_Ar,    :(species_list[sID.Ar].gamma) )

    ReplaceSymbol!(expr, :uB_Ar_m,       :(species_list[sID.Ar_m].v_Bohm) )
    ReplaceSymbol!(expr, :vth_Ar_m,      :(species_list[sID.Ar_m].v_thermal) )
    ReplaceSymbol!(expr, :D_Ar_m,        :(species_list[sID.Ar_m].D) )
    ReplaceSymbol!(expr, :gamma_Ar_m,    :(species_list[sID.Ar_m].gamma) )

    ReplaceSymbol!(expr, :uB_Ar_r,       :(species_list[sID.Ar_r].v_Bohm) )
    ReplaceSymbol!(expr, :vth_Ar_r,      :(species_list[sID.Ar_r].v_thermal) )
    ReplaceSymbol!(expr, :D_Ar_r,        :(species_list[sID.Ar_r].D) )
    ReplaceSymbol!(expr, :gamma_Ar_r,    :(species_list[sID.Ar_r].gamma) )

    ReplaceSymbol!(expr, :uB_Ar_4p,       :(species_list[sID.Ar_4p].v_Bohm) )
    ReplaceSymbol!(expr, :vth_Ar_4p,      :(species_list[sID.Ar_4p].v_thermal) )
    ReplaceSymbol!(expr, :D_Ar_4p,        :(species_list[sID.Ar_4p].D) )
    ReplaceSymbol!(expr, :gamma_Ar_4p,    :(species_list[sID.Ar_4p].gamma) )
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