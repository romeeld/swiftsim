#!/bin/bash

#SBATCH -J CW_test
#SBATCH -N 1
#SBATCH --time=12:00:00
#SBATCH --partition=cosma7
#SBATCH -A dp004
#SBATCH --exclusive

module purge
module load intel_comp/2020-update2 intel_mpi/2020-update2 ucx/1.8.1 parmetis/4.0.3-64bit parallel_hdf5/1.10.6 gsl/2.5 fftw/3.3.8cosma7
          
./swift --pin --hydro --threads=28 blob.yml

