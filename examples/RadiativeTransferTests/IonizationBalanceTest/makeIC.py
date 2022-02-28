#!/usr/bin/env python3

# -----------------------------------------------------------
# Use 8^3 particles in 3D to test the
# photo-ionization balance setup within SPHM1-RT.
# -----------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import unyt
import numpy as np

#Gadget-oid units of 10^10 Msun, Mpc, and km/s
unitsystem = cosmo_units 

# Box is 20 kpc
boxsize = 0.02 * unitsystem["length"]

# hydrogen density in cm^{-3}
nH_cgs =  1.0 / unyt.cm**3

# number of photon groups
nPhotonGroups = 4

# number of particles in each dimension
n_p = 8
# total number of particles
numPart = n_p**3

# filename of ICs to be generated
outputfilename = "ionization_balance_test.hdf5"

# index array
ids = np.linspace(1, numPart, numPart)

# particle positions
xp = unyt.unyt_array(np.zeros((n_p**3, 3), dtype=np.float64), boxsize.units)
dx = boxsize / n_p

index = 0
for i in range(n_p):
    for j in range(n_p):
        for k in range(n_p): 
            ik = int(ids[index]-1)
            xp[ik, 0] = (i + 0.5) * dx
            xp[ik, 1] = (j + 0.5) * dx 
            xp[ik, 2] = (k + 0.5) * dx  
            index += 1

# Generate object. cosmo_units corresponds to default Gadget-oid units
# of 10^10 Msun, Mpc, and km/s
w = Writer(cosmo_units, boxsize, dimension=3)

w.gas.coordinates = xp
w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
density = nH_cgs * unyt.proton_mass
volume_particle = boxsize**3 / numPart 
w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * density * volume_particle
# Generate internal energy corresponding to 10^2 K
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=np.float64) * (1e2 * unyt.kb * unyt.K) / (density * volume_particle)
)
# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
w.write(outputfilename)
