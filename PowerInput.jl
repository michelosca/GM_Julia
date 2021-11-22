module PowerInput

using SharedData: System, Species, SpeciesID

function PowerInputFunction(species::Species, system::System,
    sID::SpeciesID)

    if (species.id == sID.electron)
        S_abs = system.drivP/system.V
    else
        S_abs = 0
    end

    return S_abs
end

end