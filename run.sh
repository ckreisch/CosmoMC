#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH -n 32
#SBATCH -t 48:00:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ckreisch@astro.princeton.edu

source /home/ckreisch/.bashrc

cd /tigress/ckreisch/vanilla/CosmoMC_multinest

mpirun -n 8 --cpus-per-proc 4 ./cosmomc test.ini
