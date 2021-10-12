using Plots
using LaTeXStrings
pyplot()
p1 = plot(sol, vars=(0,2), label = L"T_e", xlabel= "Time [s]", ylabel = "Temperature [K]")
p2 = plot(sol, vars=(0,1), yscale=:log10, label = L"n_e", ylabel = latexstring("Density [m\$^{-3}\$]"))
p2 = plot!(sol, vars=(0,3), label = L"Ar")
p2 = plot!(sol, vars=(0,4), label = L"Ar^+")
p2 = plot!(sol, vars=(0,5), label = L"Ar*")
p2 = plot!(sol, vars=(0,6), label = L"Ar_2*")
p2 = plot!(sol, vars=(0,7), label = L"Ar_2^+", xlabel = "Time [s]")
plot(p1,p2, layout=(1,2))
