#!/bin/bash
#SBATCH --job-name=render_server
#SBATCH --output=render_%j.out
#SBATCH -n 16
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=15GB
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

set -e
module load cesga/2022 gcc/system openmpi/4.1.4 paraview/5.10.1

srun

