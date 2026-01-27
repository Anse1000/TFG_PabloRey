#!/bin/bash
#SBATCH -J test_cuda1
#SBATCH -t 5:00:00
#SBATCH -o test_cuda1_%j.out
#SBATCH --mem=247GB
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --gres=gpu:a100:2
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

set -e
ulimit -s unlimited
# Limpiar y cargar módulos
module purge
module load cesga/2025 nvhpc/25.3-CUDA-system

#Limpiar configuraciones previas
rm -rf build
mkdir -p build
# Compilar el programa
export HDF5_ROOT=$HOME/hdf5-nvhpc
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=ON -S . -B ./build
cmake --build ./build
# Ejecutar
./build/TFG_PRA /mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados /mnt/lustre/scratch/nlsas/home/ulc/es/pra/resultados/galactic_positions_cuda 3