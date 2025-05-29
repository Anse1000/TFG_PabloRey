#!/bin/bash
BASE_URL="https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
TABLE="I/352/gedr3dis"
FIELDS="Source,rgeo"
N_PARTS=500

# Define los límites del source_id
SOURCE_MIN=4295806720
SOURCE_MAX=6917528997577384320
STEP=$(( (SOURCE_MAX - SOURCE_MIN) / N_PARTS ))

for ((k = 0; k < N_PARTS; k++)); do
  START=$(( SOURCE_MIN + k * STEP ))
  END=$(( START + STEP - 1 ))
  FILENAME="gedr3dis_range_${k}.csv"

  if [[ -s "$FILENAME" ]]; then
    echo "Parte $k ya existe. Se omite."
    continue
  fi

  echo "Descargando parte $k: source_id entre $START y $END"

  QUERY="SELECT ${FIELDS} FROM \"${TABLE}\" WHERE Source BETWEEN ${START} AND ${END}"
  ENCODED_QUERY=$(echo "$QUERY" | sed 's/ /+/g' | sed 's/"/%22/g')
  URL="${BASE_URL}?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=${ENCODED_QUERY}"

  for attempt in {1..3}; do
    echo "Intento $attempt..."
    wget -O "$FILENAME" "$URL" && break
    echo "Falló el intento $attempt para parte $k"
    sleep 5
  done

  if [[ ! -s "$FILENAME" ]]; then
    echo "Falló la descarga de parte $k tras 3 intentos." >&2
    rm -f "$FILENAME"
  else
    echo "Parte $k descargada correctamente."
  fi
done