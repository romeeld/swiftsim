import numpy as np
import unyt


# get the unit for rt according to the RT scheme:
def get_units(scheme,unit_system="cgs_units"):
    dictunit = {}
    if unitbase == "cgs_units":
        if scheme.startswith("GEAR M1closure"):
            dictunit["energy_units"] = unyt.erg
            dictunit["energy_units_str"] = "\\rm{erg}"
            dictunit["flux_units"] = 1e10 * energy_units / unyt.cm ** 2 / unyt.s
            dictunit["flux_units_str"] = "10^{10} \\rm{erg} \\ \\rm{cm}^{-2} \\ \\rm{s}^{-1}"
        elif scheme.startswith("SPH M1closure"):
            dictunit["energy_units"] = unyt.erg
            dictunit["energy_units_str"] = "\\rm{erg}"
            dictunit["flux_units"] = 1e10 * energy_units * unyt.cm / unyt.s
            dictunit["flux_units_str"] = "10^{10} \\rm{erg} \\ \\rm{cm} \\ \\rm{s}^{-1}"        
        else:
            print("RT scheme not identified. Exit.")
            exit()
    elif unitbase == "cosmo_units":
        if scheme.startswith("GEAR M1closure"):
            dictunit["energy_units"] = unyt.erg
            dictunit["energy_units_str"] = "\\rm{erg}"
            dictunit["flux_units"] = 1e10 * energy_units / unyt.cm ** 2 / unyt.s
            dictunit["flux_units_str"] = "10^{10} \\rm{erg} \\ \\rm{cm}^{-2} \\ \\rm{s}^{-1}"
        elif scheme.startswith("SPH M1closure"):
            dictunit["energy_units"] = unyt.erg
            dictunit["energy_units_str"] = "\\rm{erg}"
            dictunit["flux_units"] = 1e10 * energy_units * unyt.cm / unyt.s
            dictunit["flux_units_str"] = "10^{10} \\rm{erg} \\ \\rm{cm} \\ \\rm{s}^{-1}"        
        else:
            print("RT scheme not identified. Exit.")
            exit()
    return dictunit
