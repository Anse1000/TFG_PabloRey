#!/bin/bash
#SBATCH -J test_cuda          # Nombre del trabajo
#SBATCH -t 10:00:00           # Tiempo máximo de ejecución (2 horas)
#SBATCH -o test_cuda.out      # Archivo de salida estándar
#SBATCH --mem=50000           # Memoria solicitada (en MB)
#SBATCH --gres=gpu:T4:1
#SBATCH --mail-type=END,TIME_LIMIT_90
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

# Cargar módulos (si es necesario)
module load gnu11/11.2.1
module load cuda/11.8.0
# Ejecutar el programa
make CUDA=1
chmod u+x TFG_PRA
./TFG_PRA pruebas