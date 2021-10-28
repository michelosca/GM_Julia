module SolveSystem

using InputBlock_Species: dens_offset
using DifferentialEquations
using Plots

function ExecuteProblem(eq_list, init, tspan)

    prob = ODEProblem(ode_fn!, init, tspan, eq_list)
    sol = solve(prob)

    p1 = plot(sol, vars=(0,1), label = "Te")
    p2 = plot(sol, vars=(0,2), label = "ne")
    display(plot(p1))
    display(plot(p2))



end

function ode_fn!(dy,y,p,t)
    # dy are the derivatives of y, i.e. the lhs of the ODE equations
    # y are [temperature,dens] (the input Vector in other functions) 
    # p are the free parameters the y-functions may have
    # t is time

    # Make sure lenght of y is the same as p
    n = length(p)

    for i in 1:n
        dy[i] = GatherListOfFunctions(y, p[i])
    end
end


function GatherListOfFunctions(input::Vector{Float64}, funct_list::Vector{Function})

    temp = input[1:dens_offset]
    dens = input[dens_offset + 1:end]

    f_out = 0
    for current_f in funct_list
        f_curr = current_f(dens, temp)
        f_out += f_curr
    end

    return f_out
end

end