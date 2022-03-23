#!/usr/bin/env python3

###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
# Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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


# -----------------------------------------------------------
# Add initial conditions for photon energies and fluxes
# for 1D diffusion of photons.
# -----------------------------------------------------------

from swiftsimio import Writer

import unyt
import numpy as np
import h5py

# define unit system to use
unitsystem = unyt.unit_systems.cgs_unit_system

# Box is 1e7 cm
boxsize = 1e7 * unitsystem["length"]

# number of photon groups
nPhotonGroups = 1

# (reduced) speed of light
cred = unyt.c

# opacity
chi_opac = 1.e-3 * unitsystem["length"]**2 / unitsystem["mass"]

# density
rho0 = 1.0 * unitsystem["mass"] / unitsystem["length"]**3 

# total radiation energy
Ecr0 = 1.0 * unitsystem["energy"]

# number of particles in each dimension
n_p = 1000

# time
time = 0.0 * unitsystem["time"] 

# width of the initial Gaussian 
offset = 0.01 * boxsize 


# filename of ICs to be generated
outputfilename = "diffusion_1D.hdf5"

# Diffusion Coefficient
# since we have not corrected for 1d in the code, we include 1/3 there
kappa_dif = cred / chi_opac / rho0 / 3.0


def diffusionsol(Ecr0,kappa,time,rneed,offset,dim=1.):
    """
    define a function to generate Gaussian (also diffusion solution):
    Ecr0              # total radiation energy
    kappa                # Diffusion Coefficient
    offset             # width of the initial Gaussian 
    """
    ecr = Ecr0/np.power(2.*np.pi*(2.*kappa*time+offset*offset),dim/2.)*np.exp(-rneed*rneed/2./(2.*kappa*time+offset*offset))
    return ecr

def ddiffusionsol(Ecr0,kappa,time,rneed,offset,dim=1.):
    decr = -rneed * Ecr0/np.power(2.*np.pi,dim/2.)/np.power(2.*kappa*time+offset*offset,dim/2.+1.0)*np.exp(-rneed*rneed/2./(2.*kappa*time+offset*offset)) 
    return  decr



def initial_condition(x, V):
    """
    The initial conditions

    x: particle position. 3D unyt array
    V: particle "volume". 1D unyt array or scalar

    returns: 
    E: photon energy density for each photon group. List of scalars with size of nPhotonGroups
    F: photon flux for each photon group. List with size of nPhotonGroups of numpy arrays of shape (3,)
    """

    # you can make the photon quantities unitless, the units will
    # already have been written down in the writer.

    E_list = []
    F_list = []

    # Group 1 Photons:
    # -------------------

    E = diffusionsol(Ecr0,kappa_dif,time,x[0]-boxsize*0.5,offset) 
    F = np.zeros(3, dtype=np.float64) * E * cred 
    # since we have not corrected for 1d in the code, we include 1/3 there
    F[0] = - cred / chi_opac / rho0 / 3.0 * ddiffusionsol(Ecr0,kappa_dif,time,x[0]-boxsize*0.5,offset)

    E_list.append(E)
    F_list.append(F)

    return E_list, F_list


if __name__ == "__main__":

    xp = unyt.unyt_array(np.zeros((n_p, 3), dtype=np.float64), boxsize.units)

    dx = boxsize / n_p

    for i in range(n_p):
        xp[i, 0] = (i + 0.5) * dx

    w = Writer(unyt.unit_systems.cgs_unit_system, boxsize, dimension=1)

    w.gas.coordinates = xp
    w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
    #we need work around for 1d case:
    m_tocode = rho0 * boxsize / n_p 
    w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * np.array(m_tocode.in_cgs()) * unyt.g
    w.gas.internal_energy = (
        np.ones(xp.shape[0], dtype=np.float64) * (300.0 * unyt.kb * unyt.K) / (unyt.g)
    )

    # Generate initial guess for smoothing lengths based on MIPS
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=1)

    # If IDs are not present, this automatically generates
    w.write(outputfilename)

    # Now open file back up again and add photon groups
    # you can make them unitless, the units have already been
    # written down in the writer. In this case, it's in cgs.

    with h5py.File(outputfilename, "r+") as F:
        header = F["Header"]
        nparts = header.attrs["NumPart_ThisFile"][0]
        parts = F["/PartType0"]

        for grp in range(nPhotonGroups):
            dsetname = "PhotonEnergiesGroup{0:d}".format(grp + 1)
            energydata = np.zeros((nparts), dtype=np.float32)
            parts.create_dataset(dsetname, data=energydata)

            dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
            #  if dsetname not in parts.keys():
            fluxdata = np.zeros((nparts, 3), dtype=np.float32)
            parts.create_dataset(dsetname, data=fluxdata)
        
        Elist = []
        Flist = []
        for p in range(nparts):
            E, Flux = initial_condition(xp[p], dx)
            for g in range(nPhotonGroups):
                Esetname = "PhotonEnergiesGroup{0:d}".format(g + 1)
                parts[Esetname][p] = E[g] * boxsize / n_p
                Fsetname = "PhotonFluxesGroup{0:d}".format(g + 1)
                parts[Fsetname][p] = Flux[g] * boxsize / n_p
                Elist.append(np.array(E[g] * boxsize / n_p))
                Flist.append(np.array(Flux[g][0] * boxsize / n_p))