#!/bin/bash

# Crear carpeta de salida si no existe
mkdir -p temp
mkdir -p ../reducidos_astro

# Función que procesa un solo archivo
procesar_archivo() {
    local file="$1"
    local index="$2"
    local clean_index=$((10#$index))
    local temp_csv="temp/temp_${clean_index}.csv"
    local reduced_file="../reducidos_astro/astro$(printf "%02d" "$clean_index").csv"

    echo "Procesando archivo: $file..."

    # Descomprimir sin borrar el .gz
    gunzip -c "$file" > "$temp_csv"

    # Procesar con awk y guardar el archivo reducido
    awk -F',' '
    BEGIN { foundHeader = 0 }
    /^#/ { next }
    !foundHeader {
        foundHeader = 1
        for (i=1; i<=NF; i++) colName[$i] = i
        print $colName["source_id"], $colName["mass_flame"], $colName["radius_flame"], $colName["logg_gspphot"]
        next
    }
    {
        print $colName["source_id"], $colName["mass_flame"], $colName["radius_flame"], $colName["logg_gspphot"]
    }' OFS=',' "$temp_csv" > "$reduced_file"

    rm -f "$temp_csv"
}

export -f procesar_archivo

# Obtener lista de archivos y numerarlos
ls *.gz | nl -v 1 -n rz | while read -r index file; do
    echo "$file $index"
done > archivos_con_indices.txt

# Ejecutar en paralelo con 4 procesos simultáneos (ajusta a tus recursos disponibles)
cat archivos_con_indices.txt | parallel -j 16 --colsep ' ' procesar_archivo {1} {2}

echo "Proceso completado. Archivos reducidos en la carpeta 'reducidos_astro/'."
