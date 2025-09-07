#!/bin/bash
#SBATCH -J test_ncu
#SBATCH -t 20:00:00
#SBATCH -o test_ncu_%j.out
#SBATCH --mem=247GB
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --gres=gpu:a100:1
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

ulimit -s unlimited
# Cargar módulos
module load cesga/2025 nvhpc/25.3-CUDA-system
# Compilar el programa
mkdir -p build
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_CUDA=ON -S . -B ./build
cmake --build ./build
# Ejecutar
ncu --set=section:summary \
    --target-processes all \
    --force-overwrite \
    --export ncu_summary \
    ./build/TFG_PRA /mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados /mnt/lustre/scratch/nlsas/home/ulc/es/pra/resultados/galactic_positions_cuda 1