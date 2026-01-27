#!/bin/bash
#SBATCH --job-name=test_pasos
#SBATCH --output=test_pasos_%j.out
#SBATCH -n 64
#SBATCH --time=8:00:00
#SBATCH --mem=180GB
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

set -e
module load hdf5

python -u test_pasos.py
