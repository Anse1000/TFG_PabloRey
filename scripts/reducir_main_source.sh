#!/bin/bash

# Crear carpeta de salida si no existe
mkdir -p reducidos

# Contador para los archivos procesados
count=1

# Procesar cada archivo .gz uno por uno
for file in *.gz; do
    echo "Procesando archivo: $file..."

    # Descomprimir (sin borrar el .gz)
    gunzip -c "$file" > temp.csv
    csv_file="temp.csv"

    # Generar nombre del archivo reducido
    reduced_file="reducidos/gaia_source$(printf "%02d" $count).csv"

    # Procesar con awk y guardar el archivo reducido
    awk -F',' '
    BEGIN { foundHeader = 0 }
    /^#/ { next }  # Ignora líneas que empiezan con #
    !foundHeader {  # Primera línea sin #
        foundHeader = 1
        for (i=1; i<=NF; i++) colName[$i] = i
        print $colName["source_id"], $colName["ra"], $colName["dec"], $colName["pmra"], \
              $colName["pmdec"], $colName["radial_velocity"], \
              $colName["phot_g_mean_mag"], $colName["bp_rp"]
        next
    }
    {
        print $colName["source_id"], $colName["ra"], $colName["dec"], $colName["pmra"], \
              $colName["pmdec"], $colName["radial_velocity"], \
              $colName["phot_g_mean_mag"], $colName["bp_rp"]
    }' OFS=',' "$csv_file" > "$reduced_file"

    # Borrar el archivo temporal CSV
    rm -f temp.csv

    # Incrementar el contador
    count=$((count + 1))
done

echo "Proceso completado. Archivos reducidos en la carpeta 'reducidos/'."