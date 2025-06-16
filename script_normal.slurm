#!/bin/bash
#SBATCH -J test_gaia
#SBATCH -t 00:50:00
#SBATCH -o test_gaia_%j.out
#SBATCH --mem=500GB
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
ulimit -s unlimited
# Cargar módulos
module load cesga/2025 intel
# Compilar el programa
mkdir -p build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=OFF -DENABLE_AVX512=OFF -S . -B ./build
cmake --build ./build
# Ejecutar
./build/TFG_PRA /mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados