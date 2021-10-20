module PowerInput

using SharedData
using SharedData: drivP, V 
using SharedData: Species

function PowerInputFunction(dens::Vector{Float64}, temp::Vector{Float64},
    species::Species)

    if (species.id == SharedData.s_electron_id)
        S_abs = drivP/V
    else
        S_abs = 0
    end

    return S_abs
end

end