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
using SharedData: p_square

function PowerInputFunction(species::Species, system::System,
    sID::SpeciesID, t_sim::Float64)

    S_abs = 0.0 
    if (species.id == sID.electron)
        if system.P_shape == p_square
            # Time period fraction
            t_start = system.P_start
            t_offset = t_sim - t_start
            freq = system.drivf
            period_fraction = t_offset * freq
            period_fraction = period_fraction - round(period_fraction)

            # Time fraction whithin the on-period
            dr = system.P_duty_ratio
            dr_ratio = period_fraction / dr
            
            # Depending at which dr_ratio there may be a power ramping up/down
            power_rate = 0.0
            if dr_ratio <= 1.0 && dr_ratio >= 0.0
                on_slope = system.on_slope
                off_slope = system.off_slope
                if dr_ratio < on_slope
                    # Ramp up at start of pulse-on
                    power_rate = dr_ratio/on_slope
                elseif dr_ratio > 1.0 - off_slope
                    # Ramp down at end of pulse-on
                    power_rate = (1.0-dr_ratio)/off_slope
                else
                    # Constant power between start and end on-pulse
                    power_rate = 1.0
                end
            end
            S_abs = system.P_absorbed * power_rate
        else
            S_abs = system.P_absorbed
        end
    end

    return S_abs
end

end