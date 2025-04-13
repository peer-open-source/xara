#!/bin/bash               
#SBATCH -J myjob                      # Job name
#SBATCH -o myjob.o%j                  # Name of stdout output file
#SBATCH -e myjob.e%j                  # Name of stderr error file
#SBATCH -p small                      # Queue (partition) name
#SBATCH -N 1                          # Total # of nodes (must be 1 for serial)
#SBATCH -n 9                          # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 20:00:00                   # Run time (hh:mm:ss)
#SBATCH --mail-type=all               # Send email at begin and end of job
#SBATCH -A DesignSafe-SimCenter       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=fmckenna@berkeley.edu

module load intel
module load petsc
module load hdf5
set -x                                 #{echo cmds, use "set echo" in csh}
ibrun OpenSeesMP fmk3.tcl              # Run the MPI executable named "a.out"
