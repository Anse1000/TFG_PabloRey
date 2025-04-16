#!/bin/bash
#SBATCH -J test_cuda
#SBATCH -t 10:00:00
#SBATCH -o test_cuda.out
#SBATCH --mem=50000
#SBATCH --gres=gpu:T4:1

# Cargar módulos (si es necesario)
module load gnu11/11.2.1
module load cuda/11.8.0
module load cmake/3.25.2
# Compilar el programa
mkdir -p build && cd build && cmake .. -DENABLE_CUDA=ON && cmake --build .
# Ejecutar el programa
./TFG_PRA ../pruebas