module SolveSystem

using DifferentialEquations: ODEProblem, solve
using DifferentialEquations: Tsit5 
using 

function ExecuteProblem(eq_list, init, tspan, offset)

    prob = ODEProblem(ode_fn!, init, tspan, (eq_list,offset))
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    return sol
end

function ode_fn!(dy,y,p,t)
    # dy are the derivatives of y, i.e. the lhs of the ODE equations
    # y are [temperature,dens] (the input Vector in other functions) 
    # p are the free parameters the y-functions may have
    # t is time

    # Make sure lenght of y is the same as p
    eq_list = p[1]
    offset = p[2]
    system = p[3]
    species_list = p[4]
    reaction_list = p[5]

    if (system.power_input_method == p_ccp_id)
        V_sheath = MeanSheathVoltage(dens, system)
    elseif (system.power_input_method == p_ccp_id)
        V_sheath = InterpolateSheathVoltage(dens, temp, species_list,
            reaction_list, system)
    end

    n = length(eq_list)
    for i in 1:n
        dy[i] = GatherListOfFunctions(y, eq_list[i], offset, V_sheath)
    end
end


function GatherListOfFunctions(input::Vector{Float64}, funct_list::Vector{Function},
    offset::Int64, V_sheath::Float64)

    temp = input[1:offset]
    dens = input[offset + 1:end]

    f_out = 0
    for current_f in funct_list
        f_curr = current_f(dens, temp, V_sheath)
        f_out += f_curr
    end

    return f_out
end

end