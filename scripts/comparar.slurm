#!/bin/bash
#SBATCH -J comparar_archivos
#SBATCH -o comparar_%j.out
#SBATCH --error=comparar_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH -t 10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

# -------------------------------
# CONFIGURACIÓN
# -------------------------------
FILE1="galactic_positions_cuda.csv"
FILE2="galactic_positions_normal.csv"

echo "Usando TMPDIR: $TMPDIR"

# -------------------------------
# 1. Ordenar los CSV por ID
# -------------------------------
csv1="$TMPDIR/csv1"
csv2="$TMPDIR/csv2"

echo "Ordenando $FILE1 ..."
# Quitamos cabecera y ordenamos por ID
( head -n1 "$FILE1" && tail -n +2 "$FILE1" | sort -T "$TMPDIR" -t, -k1,1n \
    --parallel=$SLURM_CPUS_PER_TASK -S 6G ) > "$csv1"

echo "Ordenando $FILE2 ..."
( head -n1 "$FILE2" && tail -n +2 "$FILE2" | sort -T "$TMPDIR" -t, -k1,1n \
    --parallel=$SLURM_CPUS_PER_TASK -S 6G ) > "$csv2"

# -------------------------------
# 2. Comparar archivos
# -------------------------------
echo "Comparando archivos..."
# Saltamos la cabecera y usamos join
tail -n +2 "$csv1" > "${csv1}_data"
tail -n +2 "$csv2" > "${csv2}_data"

join -t, -j 1 "${csv1}_data" "${csv2}_data" | awk -F, '{
    id=$1
    x1=$2; y1=$3; z1=$4
    x2=$5; y2=$6; z2=$7

    dx=x2-x1
    dy=y2-y1
    dz=z2-z1

    if (dx!=0 || dy!=0 || dz!=0) {
        printf "ID:%s ΔX=%s ΔY=%s ΔZ=%s\n", id, dx, dy, dz
    }
}' > diferencias.txt

# -------------------------------
# 3. IDs únicos
# -------------------------------
echo "Buscando IDs únicos..."
cut -d, -f1 "${csv1}_data" > ids1.txt
cut -d, -f1 "${csv2}_data" > ids2.txt

comm -23 ids1.txt ids2.txt > solo_en_archivo1.txt
comm -13 ids1.txt ids2.txt > solo_en_archivo2.txt

# -------------------------------
# 4. Reporte final
# -------------------------------
echo "Proceso terminado:"
echo " - diferencias.txt"
echo " - solo_en_archivo1.txt"
echo " - solo_en_archivo2.txt"

