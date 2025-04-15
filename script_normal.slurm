#!/bin/bash
#SBATCH -J test_gaia
#SBATCH -t 70:00:00
#SBATCH -o test_gaia.out
#SBATCH --mem=50000
#SBATCH --partition=compute3

# Cargar módulos (si es necesario)
module load gnu11/11.2.1
# Ejecutar el programa
make
chmod u+x TFG_PRA
./TFG_PRA pruebas