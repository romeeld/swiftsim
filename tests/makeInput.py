###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

import h5py
from numpy import *

# Generates a swift IC file containing a cartesian distribution of particles
# at a constant density and pressure in a cubic box

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1.0
L = 4  # Number of particles along one axis
density = 2.0  # Density
P = 1.0  # Pressure
gamma = 5.0 / 3.0  # Gas adiabatic index
material = 0  # Ideal gas
fileName = "input.hdf5"

# ---------------------------------------------------
numPart = L ** 3
mass = boxSize ** 3 * density / numPart
internalEnergy = P / ((gamma - 1.0) * density)

# chemistry data
he_density = density * 0.24

# Generate particles
coords = zeros((numPart, 3))
v = zeros((numPart, 3))
m = zeros((numPart, 1))
h = zeros((numPart, 1))
u = zeros((numPart, 1))
rho = zeros((numPart, 1))
ids = zeros((numPart, 1), dtype="L")
mat = zeros((numPart, 1), dtype="i")

# chemistry data
he = zeros((numPart, 1))

for i in range(L):
    for j in range(L):
        for k in range(L):
            index = i * L * L + j * L + k
            x = i * boxSize / L + boxSize / (2 * L)
            y = j * boxSize / L + boxSize / (2 * L)
            z = k * boxSize / L + boxSize / (2 * L)
            coords[index, 0] = x
            coords[index, 1] = y
            coords[index, 2] = z
            v[index, 0] = 0.0
            v[index, 1] = 0.0
            v[index, 2] = 0.0
            m[index] = mass
            h[index] = 2.251 * boxSize / L
            u[index] = internalEnergy
            rho[index] = density
            ids[index] = index
            mat[index] = material
            # chemistry data
            he[index] = he_density


# --------------------------------------------------

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]


# Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0


# Particle group
grp = file.create_group("/PartType0")
ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
ds[()] = coords
ds = grp.create_dataset("Velocities", (numPart, 3), "f")
ds[()] = v
ds = grp.create_dataset("Masses", (numPart, 1), "f")
ds[()] = m
ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
ds[()] = h
ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
ds[()] = u
ds = grp.create_dataset("Density", (numPart, 1), "f")
ds[()] = rho
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids
ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
ds[()] = mat
# chemistry
ds = grp.create_dataset("HeDensity", (numPart, 1), "f")
ds[()] = he

file.close()
