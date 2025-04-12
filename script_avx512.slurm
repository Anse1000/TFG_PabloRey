#!/bin/bash
#SBATCH -J test_gaia
#SBATCH -t 20:00:00
#SBATCH -o test_gaia.out
#SBATCH --mem=50000
#SBATCH --mail-type=END,TIME_LIMIT_90
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

# Cargar módulos (si es necesario)
module load gnu11/11.2.1
# Ejecutar el programa
make AVX_512=1
chmod u+x TFG_PRA
./TFG_PRA pruebas