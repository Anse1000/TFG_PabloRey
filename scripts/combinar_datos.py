import pandas as pd
import glob
import os
from concurrent.futures import ProcessPoolExecutor

# Carpetas de entrada/salida
astro_folder = 'astro'
mainsource_folder = 'reducidos'
distance_folder = 'distances'
output_folder = 'optimizados'

# Asegurarse de que existe la carpeta de salida
os.makedirs(output_folder, exist_ok=True)

# Listar y ordenar archivos
mainsource_files = sorted(glob.glob(os.path.join(mainsource_folder, '*.csv')))
astro_files = sorted(glob.glob(os.path.join(astro_folder, '*.csv')))
distance_files = sorted(glob.glob(os.path.join(distance_folder, '*.csv')))

# Cargar todos los archivos de distancia en uno solo
distance_df = pd.concat([pd.read_csv(f) for f in distance_files], ignore_index=True)
# Renombrar columna 'Source' a 'source_id' si existe
if 'Source' in distance_df.columns:
    distance_df.rename(columns={'Source': 'source_id'}, inplace=True)

# Función de mezcla individual
def merge_and_save(i, main_file, astro_file, distance_df, output_folder):
    astro_df = pd.read_csv(astro_file)
    main_df = pd.read_csv(main_file)

    # Eliminar columna 'parallax' si existe
    if 'parallax' in main_df.columns:
        main_df.drop(columns=['parallax'], inplace=True)

    merged = pd.merge(main_df, astro_df, on='source_id', how='inner')
    merged = pd.merge(merged, distance_df, on='source_id', how='left')

    output_file = os.path.join(output_folder, f'merged_{i:04d}.csv')
    merged.to_csv(output_file, index=False)
    return f'Guardado: {output_file}'

# Empaquetar argumentos para ejecución paralela
def run_parallel_merge():
    with ProcessPoolExecutor(max_workers=16) as executor:
        tasks = [
            executor.submit(merge_and_save, i, main_file, astro_file, distance_df, output_folder)
            for i, (main_file, astro_file) in enumerate(zip(mainsource_files, astro_files))
        ]

        for future in tasks:
            print(future.result())


if __name__ == '__main__':
    run_parallel_merge()
