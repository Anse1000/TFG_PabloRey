#!/bin/bash
#SBATCH -J comparar_archivos
#SBATCH -o comparar_%j.out
#SBATCH --error=comparar_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=128G               # Ajusta según RAM disponible
#SBATCH -t 10:00:00          # Tiempo máximo (48h)
#SBATCH --mail-type=END,FAIL     # Opcional: avisos por correo
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

# -------------------------------
# CONFIGURACIÓN
# -------------------------------
FILE1="galactic_positions_cuda.txt"
FILE2="galactic_positions_normal.txt"

echo "Usando TMPDIR: $TMPDIR"

# -------------------------------
# 1. Convertir y ordenar externamente
# -------------------------------
csv1="$TMPDIR/csv1"
csv2="$TMPDIR/csv2"

echo "Convirtiendo y ordenando $FILE1 ..."
sed -E 's/ID: ([0-9]+) X: ([^ ]+) Y = ([^,]+), Z = (.+)/\1,\2,\3,\4/' "$FILE1" \
| sort -T "$TMPDIR" -t, -k1,1 --parallel=$SLURM_CPUS_PER_TASK -S 6G > "$csv1"

echo "Convirtiendo y ordenando $FILE2 ..."
sed -E 's/ID: ([0-9]+) X: ([^ ]+) Y = ([^,]+), Z = (.+)/\1,\2,\3,\4/' "$FILE2" \
| sort -T "$TMPDIR" -t, -k1,1 --parallel=$SLURM_CPUS_PER_TASK -S 6G > "$csv2"

# -------------------------------
# 2. Comparar archivos
# -------------------------------
echo "Comparando archivos..."
join -t, -j 1 "$csv1" "$csv2" | awk -F, '{
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
comm -23 <(cut -d, -f1 "$csv1") <(cut -d, -f1 "$csv2") > solo_en_archivo1.txt
comm -13 <(cut -d, -f1 "$csv1") <(cut -d, -f1 "$csv2") > solo_en_archivo2.txt

# -------------------------------
# 4. Reporte final
# -------------------------------
echo "Proceso terminado:"
echo " - diferencias.txt"
echo " - solo_en_archivo1.txt"
echo " - solo_en_archivo2.txt"
