.. SPHM1RT Radiative Transfer
    Tsang Keung Chan 01.2022

.. _rt_SPHM1:
   
SPHM1 RT
-------

.. warning::
    The radiative transfer schemes are still in development and are not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.


Compiling for SPHM1-RT
~~~~~~~~~~~~~~~~~~~~~

-   To compile swift to be able to run with SPHM1-RT, you need to configure with
    ``--with-rt=SPHM1RT_N`` where ``N`` is the integer number of photon groups that 
    you intend to use in your simulation.

-   You need to choose a Riemann solver for the RT equations. You can choose
    between the ``GLF`` and ``HLL`` solver. For the time being, I recommend 
    sticking to the ``GLF`` solver as the ``HLL`` solver is more expensive,
    but seemingly offers no advantage, although this remains to be comfirmed
    in further testing.

-   SPHM1-RT is compatible with any SPH scheme. You'll
    need to compile using ``--with-hydro=sphenix`` or other SPH schemes, e.g. we have tested gadget2, minimal, and sphenix.




Compulsory Runtime Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to provide the following runtime parameters in the yaml file:

.. code:: yaml

SPHM1RT:
    cred: 2.99792458e10                                 # value of reduced speed of light for the RT solver in code unit
    CFL_condition: 0.1                                  # CFL condition for RT, independent of hydro 
    chi:  [0, 0, 0]                                     # (Optional) initial opacity in code unit for all gas particles
    photon_groups_Hz: [3.288e15, 5.945e15, 13.157e15]   # Photon frequency group bin edges in Hz. Needs to be 1 less than the number of groups (N) requested during the configuration (--with-RT=SPHM1RT_N). Outer edges of zero and infinity are assumed.


The ``photon_groups_Hz`` need to be ``N - 1`` frequency edges (floats) to separate 
the spectrum into ``N`` groups. The outer limits of zero and infinity are 
assumed.

Initial Conditions
~~~~~~~~~~~~~~~~~~


Setting Up Initial Conditions for RT
````````````````````````````````````

Optionally, you may want to provide initial conditions for the radiation field
and/or the mass fraction of the ionizing species.
To do so, you need to add the following datasets to the ``/PartType0`` particle
group:

.. code:: 

   PhotonEnergiesGroup1
   PhotonEnergiesGroup2 
   .
   .
   .
   PhotonEnergiesGroupN
   PhotonFluxesGroup1
   PhotonFluxesGroup2
   .
   .
   .
   PhotonFluxesGroupN


The ``PhotonEnergies*`` datasets need to have dimension ``nparts``, while the
``PhotonFluxesGroup*`` datasets need to have dimension ``(nparts, 3)``, where
``nparts`` is the number of hydro particles. If you are writing initial
conditions where the fields have units, then ``PhotonEnergies*`` are expected to
have units of energy :math:`[M L^2 T^{-2}]`), while the ``PhotonFluxes*`` fields
should be in units of energy flux (energy times speed), :math:`[M L^3
T^{-3}]`).


Example using Python and ``swiftsimio``
````````````````````````````````````````

If you are using `swiftsimio <https://github.com/SWIFTSIM/swiftsimio>`_ to write
the initial condition files, then the easiest way of adding the RT initial
conditions is to first use the swiftsimio routines to write a file, then open it
up again and write the additional RT fields again using ``h5py`` routines.

Here is an example:

.. code:: python

    from swiftsimio import Writer
    import unyt
    import numpy as np
    import h5py

    # define unit system to use.
    unitsystem = unyt.unit_systems.cgs_unit_system

    # number of photon groups
    nPhotonGroups = 4

    # filename of ICs to be generated
    outputfilename = "my_rt_ICs.hdf5"

    # open a swiftsimio.Writer object
    w = Writer(...)

    # do your IC setup for gas, gravity etc now
    # ... 

    # write the IC file without doing anything RT related.
    w.write(outputfilename)

    # Now open file back up again and add RT data.
    F = h5py.File(outputfilename, "r+")
    header = F["Header"]
    nparts = header.attrs["NumPart_ThisFile"][0]
    parts = F["/PartType0"]

    # Create initial photon energies and fluxes. You can leave them unitless, 
    # the units have already been written down with w.write(). In this case, 
    # it's in cgs.
    for grp in range(nPhotonGroups):
        dsetname = "PhotonEnergiesGroup{0:d}".format(grp + 1)
        energydata = np.ones((nparts), dtype=np.float32) * some_value_you_want
        parts.create_dataset(dsetname, data=energydata)

        dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
        fluxdata = np.zeros((nparts, 3), dtype=np.float32) * some_value_you_want
        parts.create_dataset(dsetname, data=fluxdata)

    # close up, and we're done!
    F.close()
