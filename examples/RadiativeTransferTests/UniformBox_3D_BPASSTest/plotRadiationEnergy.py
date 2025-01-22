#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

# ----------------------------------------------------
# plots 2D projection of radiation energy and fluxes
#
# Usage: give snapshot number as cmdline arg to plot
# single snapshot, otherwise this script plots
# all snapshots available in the workdir.
# ----------------------------------------------------

import gc
import os
import sys
import warnings

import matplotlib as mpl
import numpy as np
import swiftsimio
import unyt
import h5py
from matplotlib import pyplot as plt
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Parameters users should/may tweak

# plot all groups and all photon quantities
plot_all_data = True
# snapshot basename
snapshot_base = "output"
# fancy up the plots a bit?
fancy = True

# parameters for imshow plots
imshow_kwargs = {"origin": "lower", "cmap": "viridis"}

# parameters for swiftsimio projections
projection_kwargs = {"resolution": 1024, "parallel": True}

# arguments for plot of references
referenceplotkwargs = {"color": "grey", "lw": 4, "alpha": 0.6}

# arguments for legends
legendprops = {"size": 8}

# Set Units of your choice
energy_units = unyt.erg
energy_units_str = "\\rm{erg}"
flux_units = 1e6 * energy_units / unyt.cm ** 2 / unyt.s
flux_units_str = "10^{6} \\rm{erg} \\ \\rm{cm}^{-2} \\ \\rm{s}^{-1}"
time_units = unyt.s

# Set photon escape fraction
f_esc = float(sys.argv[1])

# Configure warnings to be treated as errors
warnings.simplefilter("error")

# Check the value is valid or not
if (f_esc < 0 or f_esc > 1):
    warnings.warn("f_esc is invalid; check your logic!", UserWarning)

# -----------------------------------------------------------------------


# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[2])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted
    and return their names as list
    """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_photons(filename, energy_boundaries=None, flux_boundaries=None):
    """
    Create the actual plot.

    filename: file to work with
    energy_boundaries:  list of [E_min, E_max] for each photon group.
                        If none, limits are set automatically.
    flux_boundaries:    list of [F_min, F_max] for each photon group.
                        If none, limits are set automatically.
    """

    print("working on", filename)

    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata

    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"][0])
    xlabel_units_str = meta.boxsize.units.latex_representation()

    #average_photon_energy=[2.97254e-11, 6.01548e-11, 2.83711e-10] * unyt.erg
    average_photon_energy=[2.864693e-11, 4.955534e-11, 8.834406e-11] * unyt.erg

    n_photon=np.zeros(3)


    for g in range(ngroups):
        # workaround to access named columns data with swiftsimio visualisaiton
        # add mass weights to remove surface density dependence in images
        new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
        en = np.float128(getattr(data.gas.photon_energies, "group" + str(g + 1)))

        en.convert_to_units(energy_units)
          
        sum_en = np.sum(en)

        n_photon[g] = sum_en / average_photon_energy[g]

        #print('photon groups:', g)
        #print(en)
        #print(sum_en)
        #print(n_photon)

    return meta.time.to("Myr"), n_photon


def plot_emitted_photon(snaplist, energy_boundaries=None, flux_boundaries=None):
    
    timelist = []

    photongroup1 = []
    photongroup2 = []
    photongroup3 = []
    #get the bin of the n_photon
    for f in snaplist:
        time, n_photons = plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
        
        timelist.append(time)
        photongroup1.append(n_photons[0])
        photongroup2.append(n_photons[1])
        photongroup3.append(n_photons[2])

    # Open the HDF5 file
    with h5py.File('BPASS_chab100.h5', 'r') as f:
        HI = np.array(f['Table_HI/block0_values'])
        HeI = np.array(f['Table_HeI/block0_values'])
        HeII = np.array(f['Table_HeII/block0_values'])

    age_100Myr = 10**(np.arange(0, 2.1, 0.1))
    age_100Myr_with_zero = np.insert(age_100Myr, 0, 0)

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)

    ax1.semilogy(age_100Myr_with_zero, f_esc * HI[3,:], label="reference", **referenceplotkwargs)
    ax1.hlines(y=f_esc * HI[3,-1], xmin=100, xmax=timelist[-1], **referenceplotkwargs)
    ax1.semilogy(timelist,photongroup1, label="obtained results")
    ax1.set_ylim(1e55, 1e62)
    ax1.set_xlabel("time [Myr]")
    ax1.set_ylabel("photon number")
    ax1.set_title("Photon group 1")
    ax1.legend(prop=legendprops)
    ax1.grid()

    ax2.semilogy(age_100Myr_with_zero, f_esc * HeI[3,:], label="reference", **referenceplotkwargs)
    ax2.hlines(y=f_esc * HeI[3,-1], xmin=100, xmax=timelist[-1], **referenceplotkwargs)
    ax2.semilogy(timelist,photongroup2, label="obtained results")
    ax2.set_ylim(1e55, 1e62)
    ax2.set_xlabel("time [Myr]")
    ax2.set_ylabel("photon number")
    ax2.set_title("Photon group 2")
    ax2.legend(prop=legendprops)
    ax2.grid()

    ax3.semilogy(age_100Myr_with_zero, f_esc * HeII[3,:], label="reference", **referenceplotkwargs)
    ax3.hlines(y=f_esc * HeII[3,-1], xmin=100, xmax=timelist[-1], **referenceplotkwargs)
    ax3.semilogy(timelist,photongroup3, label="obtained results")
    ax3.set_ylim(1e55, 1e62)
    ax3.set_xlabel("time [Myr]")
    ax3.set_ylabel("photon number")
    ax3.set_title("Photon group 3")
    ax3.legend(prop=legendprops)
    ax3.grid()




    plt.tight_layout()
    plt.savefig("Bpass_emitted_photon.png")

    


      


if __name__ == "__main__":
    
    snaplist = get_snapshot_list(snapshot_base)

    energy_boundaries = None
    flux_boundaries = None

    plot_emitted_photon(snaplist, energy_boundaries, flux_boundaries)
