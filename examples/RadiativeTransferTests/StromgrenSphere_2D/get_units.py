import numpy as np
import unyt


# get the unit for rt according to the RT scheme:
def get_units(scheme, unit_system="cgs_units"):
    if unit_system == "cgs_units":
        time_units = unyt.s
        energy_units = unyt.erg
        energy_units_str = "\\rm{erg}"
        if scheme.startswith("GEAR M1closure"):
            flux_units = 1e10 * energy_units / unyt.cm ** 2 / unyt.s
            flux_units_str = "10^{10} \\rm{erg} \\ \\rm{cm}^{-2} \\ \\rm{s}^{-1}"
        elif scheme.startswith("SPH M1closure"):
            flux_units = 1e10 * energy_units * unyt.cm / unyt.s
            flux_units_str = "10^{10} \\rm{erg} \\ \\rm{cm} \\ \\rm{s}^{-1}"
        else:
            print("RT scheme not identified. Exit.")
            exit()
    elif unit_system == "stromgren_units":
        time_units = unyt.Myr
        energy_units = 1e50 * unyt.erg
        energy_units_str = "10^{50} \\rm{erg}"
        if scheme.startswith("GEAR M1closure"):
            flux_units = 1e50 * unyt.erg / unyt.kpc ** 2 / unyt.Gyr
            flux_units_str = "10^{60} \\rm{erg} \\ \\rm{kpc}^{-2} \\ \\rm{Gyr}^{-1}"
        elif scheme.startswith("SPH M1closure"):
            flux_units = 1e50 * unyt.erg * unyt.kpc / unyt.Gyr
            flux_units_str = "10^{60} \\rm{erg} \\ \\rm{kpc} \\ \\rm{Gyr}^{-1}"
        else:
            print("RT scheme not identified. Exit.")
            exit()
    else:
        print("Unit system not identified. Exit.")
        exit()
    return time_units, energy_units, energy_units_str, flux_units, flux_units_str
