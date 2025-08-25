#!/bin/bash
#SBATCH -J test_cuda
#SBATCH -t 02:00:00
#SBATCH -o test_cuda_%j.out
#SBATCH --mem=400GB
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --gres=gpu:a100:8
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

ulimit -s unlimited
# Cargar módulos
module load cesga/2025 nvhpc/25.3-CUDA-system
# Compilar el programa
mkdir -p build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=ON -S . -B ./build
cmake --build ./build
# Ejecutar
./build/TFG_PRA /mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados /mnt/netapp2/Store_uni/home/ulc/es/pra/initial_positions_cuda /mnt/netapp2/Store_uni/home/ulc/es/pra/final_positions_cuda 1