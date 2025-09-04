#!/bin/bash
#SBATCH -J comp_csv
#SBATCH -o resultados_compare/comp_%j.out
#SBATCH --array=0-24          # Número de chunks
#SBATCH -t 01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G

module load cesga/2022 gcccore/system pandas/2.2.1

# Definir paths
ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
INPUT_DIR1="/mnt/netapp2/Store_uni/home/ulc/es/pra/galactic_positions_normal/step_1/cpu_results_${ID}.csv"
INPUT_DIR2="/mnt/netapp2/Store_uni/home/ulc/es/pra/galactic_positions_cuda/step_1/cuda_results_${ID}.csv"
OUTPUT_DIR="resultados_compare/resultados_${ID}.csv"

mkdir -p resultados_compare

srun python compare.py "$INPUT_DIR1" "$INPUT_DIR2" "$OUTPUT_DIR"



