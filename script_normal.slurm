#!/bin/bash
#SBATCH -J test_gaia
#SBATCH -t 4-00:00:00
#SBATCH -o test_gaia_%j.out
#SBATCH --mem=247GB
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=pablo.rey.ansemil@udc.es
ulimit -s unlimited
# Cargar módulos
module load cesga/2025 intel
# Compilar el programa
mkdir -p build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=OFF -S . -B ./build
cmake --build ./build
# Ejecutar
./build/TFG_PRA /mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados /mnt/netapp2/Store_uni/home/ulc/es/pra/initial_positions_normal.txt /mnt/netapp2/Store_uni/home/ulc/es/pra/galactic_positions_normal.txt 1