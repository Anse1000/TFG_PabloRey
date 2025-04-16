#!/bin/bash
SERVER="pluton.dec.udc.es"
USER="pablo.rey"

FILES=("src" "CMakeLists.txt" "script_cuda.slurm" "script_normal.slurm" "script_avx512.slurm")  # Archivos y carpetas a subir

REMOTE_DIR="STORE/TFG_PabloRey"         # Cambia esto a la ruta de destino en el servidor

echo "Conectando a $USER@$SERVER y subiendo archivos..."

{
  echo "cd $REMOTE_DIR"
  for file in "${FILES[@]}"; do
    if [ -d "$file" ]; then
      echo "put -r $file"
    else
      echo "put $file"
    fi
  done
  echo "bye"
} > sftp_batch.txt

sftp -b sftp_batch.txt $USER@$SERVER > /dev/null 2>&1

rm sftp_batch.txt

echo "Archivos subidos con éxito."
