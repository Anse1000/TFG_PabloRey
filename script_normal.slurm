#!/bin/bash
#SBATCH -J test_gaia
#SBATCH -t 70:00:00
#SBATCH -o test_gaia.out
#SBATCH --mem=250000
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=compute3

# Cargar módulos (si es necesario)
module load gnu11/11.2.1
module load cmake/3.25.2
# Compilar el programa
mkdir -p build && cd build && cmake .. && cmake --build .
./TFG_PRA ../reducidos