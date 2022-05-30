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

module PowerInput

using SharedData: System, Species, SpeciesID

function PowerInputFunction(species::Species, system::System,
    sID::SpeciesID, t_sim::Float64)

    if (species.id == sID.electron)
        if (system.P_shape == "sinusoidal")
            S_abs = system.drivP/system.V
        elseif (system.P_shape == "square")
                t_frac = t_sim * system.drivf - floor(t_sim * system.drivf)
                if (t_frac <= system.P_duty_ratio)
                    S_abs = system.drivP / system.V
                else
                    S_abs = 0.0
                end
        end
    else
        S_abs = 0.0
    end

    return S_abs
end

end