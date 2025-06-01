import pandas as pd
import duckdb
import glob
import os
from concurrent.futures import ProcessPoolExecutor
import time

# Configuración de paths
astro_folder = 'reducidos_astro'
mainsource_folder = 'reducidos'
distance_folder = 'distances'
output_folder = 'optimizados'
db_path = 'distance_data.duckdb'

# Crear carpeta de salida
os.makedirs(output_folder, exist_ok=True)

def crear_base_duckdb(distance_folder, db_path):
    print("[INFO] Creando base DuckDB desde archivos de distancia...")
    distance_files = sorted(glob.glob(os.path.join(distance_folder, '*.csv')))
    if not distance_files:
        raise FileNotFoundError("[ERROR] No se encontraron archivos en 'distances'.")

    chunks = []
    for i, file in enumerate(distance_files):
        print(f"[INFO] Leyendo archivo {i+1}/{len(distance_files)}: {file}")
        df = pd.read_csv(file)
        chunks.append(df)

    distance_df = pd.concat(chunks, ignore_index=True)
    if 'Source' in distance_df.columns:
        print("[INFO] Renombrando columna 'Source' a 'source_id'")
        distance_df.rename(columns={'Source': 'source_id'}, inplace=True)

    con = duckdb.connect(db_path)
    con.execute("DROP TABLE IF EXISTS distance")
    con.register("df", distance_df)
    con.execute("CREATE TABLE distance AS SELECT * FROM df")
    con.close()
    print("[OK] Base de datos DuckDB creada")

def merge_and_save(i, main_file, astro_file, db_path, output_folder):
    try:
        print(f"[{i:04d}] Procesando archivos: {os.path.basename(main_file)}, {os.path.basename(astro_file)}")
        main_df = pd.read_csv(main_file)
        astro_df = pd.read_csv(astro_file)

        if 'parallax' in main_df.columns:
            main_df.drop(columns=['parallax'], inplace=True)

        merged = pd.merge(main_df, astro_df, on='source_id', how='inner')

        source_ids = tuple(merged['source_id'].unique())
        if len(source_ids) == 1:
            source_ids = f"({source_ids[0]})"
        else:
            source_ids = str(source_ids)

        con = duckdb.connect(db_path, read_only=True)
        distance_df = con.execute(
            f"SELECT * FROM distance WHERE source_id IN {source_ids}"
        ).fetchdf()
        con.close()

        merged = pd.merge(merged, distance_df, on='source_id', how='left')

        output_file = os.path.join(output_folder, f'merged_{i:04d}.csv')
        merged.to_csv(output_file, index=False)
        return f'[OK] Guardado: {output_file}'

    except Exception as e:
        return f'[ERROR] {i:04d}: {str(e)}'

def run():
    start_time = time.time()
    print("[INFO] Iniciando procesamiento completo...")

    # Paso 1: Crear DB
    crear_base_duckdb(distance_folder, db_path)

    # Paso 2: Preparar listas de archivos
    main_files = sorted(glob.glob(os.path.join(mainsource_folder, '*.csv')))
    astro_files = sorted(glob.glob(os.path.join(astro_folder, '*.csv')))

    if len(main_files) != len(astro_files):
        raise ValueError("[ERROR] Diferente número de archivos en main y astro")

    print(f"[INFO] Total de pares a procesar: {len(main_files)}")

    # Paso 3: Merge en paralelo
    with ProcessPoolExecutor(max_workers=16) as executor:
        tasks = [
            executor.submit(merge_and_save, i, main_file, astro_file, db_path, output_folder)
            for i, (main_file, astro_file) in enumerate(zip(main_files, astro_files))
        ]
        for future in tasks:
            print(future.result())

    # Paso 4: Eliminar base de datos temporal
    if os.path.exists(db_path):
        os.remove(db_path)
        print(f"[OK] Base de datos eliminada: {db_path}")

    print(f"[OK] Proceso finalizado en {time.time() - start_time:.2f} segundos")

if __name__ == '__main__':
    run()
