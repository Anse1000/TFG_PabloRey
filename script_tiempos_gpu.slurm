#!/bin/bash
#SBATCH -J test_cuda
#SBATCH -t 20:00:00
#SBATCH -o test_cuda_%j.out
#SBATCH --mem=247GB
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

# Rutas
INPUT_DIR="/mnt/lustre/scratch/nlsas/home/ulc/es/pra/optimizados"
OUTPUT_DIR="/mnt/lustre/scratch/nlsas/home/ulc/es/pra/resultados/galactic_positions_cuda"

# Archivo donde se guardarán los tiempos
TIME_LOG="tiempos_gpus.txt"
echo "GPUs, Tiempo (HH:MM:SS)" > $TIME_LOG

echo "=== INICIO EXPERIMENTO ==="
date
echo ""

for NGPU in {1..8}; do
    DEVICES=$(seq -s "," 0 $((NGPU-1)))
    export CUDA_VISIBLE_DEVICES=$DEVICES

    echo "Ejecutando con $NGPU GPU(s)..."

    START=$(date +%s)
    ./build/TFG_PRA "$INPUT_DIR" "$OUTPUT_DIR" 1
    END=$(date +%s)

    ELAPSED=$((END-START))

    # Convertir a formato hh:mm:ss
    HH=$(printf "%02d" $((ELAPSED/3600)))
    MM=$(printf "%02d" $(((ELAPSED%3600)/60)))
    SS=$(printf "%02d" $((ELAPSED%60)))
    TIME_FORMAT="$HH:$MM:$SS"

    echo "$NGPU, $TIME_FORMAT" >> $TIME_LOG
    echo "Tiempo: $TIME_FORMAT"
    echo ""
done

echo "=== FIN EXPERIMENTO ==="
date
echo "Resultados guardados en: $TIME_LOG"
