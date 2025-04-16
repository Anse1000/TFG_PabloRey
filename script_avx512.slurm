#!/bin/bash
#SBATCH -J test_avx
#SBATCH -t 20:00:00
#SBATCH -o test_avx.out
#SBATCH --mem=50000
#SBATCH --partition=compute3

# Cargar módulos
module load gnu11/11.2.1
module load cmake/3.25.2
# Compilar el programa
mkdir -p build && cd build && cmake .. -DENABLE_AVX512=ON && cmake --build .
./TFG_PRA ../pruebas