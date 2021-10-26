module PowerInput

using SharedData: System, Species
using InputBlock_Species: s_electron_id

function PowerInputFunction(dens::Vector{Float64}, temp::Vector{Float64},
    species::Species, system::System)

    if (species.id == s_electron_id)
        S_abs = system.drivP/system.V
    else
        S_abs = 0
    end

    return S_abs
end

end