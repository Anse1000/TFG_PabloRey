#!/bin/bash
#SBATCH -J test_gaia
#SBATCH -t 00:15:00
#SBATCH -o test_gaia_%j.out
#SBATCH --mem=170000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64

# Cargar módulos
module load cesga/2025 intel
# Compilar el programa
mkdir -p build
cmake -DCMAKE_BUILD_TYPE=Release -S . -B ./build
cmake --build ./build
# Ejecutar
./build/TFG_PRA /mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados