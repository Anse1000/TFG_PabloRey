import h5py
import numpy as np
import xml.etree.ElementTree as ET
import os
import time
import glob
from concurrent.futures import ThreadPoolExecutor

# --- CONFIGURACION ---
# Ruta base donde se encuentran las carpetas step_0000, step_0001, etc.
BASE_DIR = "/mnt/lustre/scratch/nlsas/home/ulc/es/pra/resultados/galactic_positions_cuda/"
MAX_WORKERS = 32
# Umbral de movimiento para detectar saltos anomalos
MAX_ALLOWED_DISPLACEMENT = 10.0

def get_chunk_info(xdmf_path):
    """Parsea el archivo XDMF para identificar los fragmentos H5"""
    try:
        tree = ET.parse(xdmf_path)
        root = tree.getroot()
        base_dir = os.path.dirname(xdmf_path)
        grids = root.findall(".//Grid[@GridType='Uniform']")
        chunks = []
        total_elements = 0
        for grid in grids:
            mapping = {}
            n_elements = int(grid.find("Topology").get("NumberOfElements"))
            for attr in grid.findall("Attribute") + [grid.find("Geometry")]:
                di = attr.find("DataItem")
                if di is not None:
                    name = attr.get('Name') or 'XYZ'
                    h5_file, ds_path = di.text.strip().split(':')
                    mapping[name] = (os.path.join(base_dir, h5_file), ds_path)
            chunks.append({'mapping': mapping, 'n': n_elements, 'start': total_elements})
            total_elements += n_elements
        return chunks, total_elements
    except Exception as e:
        print(f"Error critico leyendo metadata de {xdmf_path}: {e}")
        return None, 0

def read_h5_chunk(chunk, out_ids, out_pos):
    """Lee un fragmento especifico de un archivo H5"""
    m = chunk['mapping']
    start, end = chunk['start'], chunk['start'] + chunk['n']
    try:
        with h5py.File(m['ID'][0], 'r') as f:
            out_ids[start:end] = f[m['ID'][1]][:]
            out_pos[start:end] = f[m['XYZ'][1]][:]
    except Exception as e:
        print(f"Error leyendo fragmento H5 {m['ID'][0]}: {e}")

def load_step(step_path):
    """Carga y ordena los datos de un paso temporal completo"""
    t_start_load = time.time()
    print(f"Analizando: {os.path.basename(step_path)}", flush=True)

    chunks, total_n = get_chunk_info(step_path)
    if chunks is None or total_n == 0:
        print(f"Carga abortada: No se encontraron datos.")
        return None

    print(f"  - Total particulas detectadas: {total_n:,}")
    print(f"  - Fragmentos H5 a procesar: {len(chunks)}")

    ids = np.empty(total_n, dtype=np.uint64)
    pos = np.empty((total_n, 3), dtype=np.float64)

    # Lectura paralela
    print(f"  - Iniciando lectura paralela ({MAX_WORKERS} hilos)...", flush=True)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for chunk in chunks:
            executor.submit(read_h5_chunk, chunk, ids, pos)

    t_io_done = time.time()
    print(f"  - Lectura E/S completada en {(t_io_done - t_start_load):.2f}s", flush=True)

    # Ordenacion para alinear particulas por ID
    print(f"  - Iniciando reordenamiento de datos por ID...", flush=True)
    t_sort_start = time.time()
    idx = np.argsort(ids)

    data_sorted = {
        'ids': ids[idx],
        'pos': pos[idx]
    }

    t_done = time.time()
    print(f"  - Ordenacion completada en {(t_done - t_sort_start):.2f}s", flush=True)
    print(f"  - Tiempo total de carga del paso: {(t_done - t_start_load):.2f}s", flush=True)

    return data_sorted

def analyze_motion():
    # Obtener lista de archivos ordenados
    step_files = sorted(glob.glob(os.path.join(BASE_DIR, "step_*/step_*.xmf")))

    if len(step_files) < 2:
        print(f"Error: Se encontraron {len(step_files)} archivos. Se requieren al menos 2 para comparar movimiento.")
        return

    print("="*70)
    print(f"INICIANDO ANALISIS TEMPORAL DE {len(step_files)} PASOS")
    print("="*70)

    # Cargar primer paso para comparar
    prev_step_data = load_step(step_files[0])
    if prev_step_data is None: return

    print("-" * 50)

    for i in range(1, len(step_files)):
        current_step_path = step_files[i]
        current_step_data = load_step(current_step_path)

        if current_step_data is None:
            print(f"Saltando comparacion para {current_step_path} debido a errores de carga.")
            continue

        # Verificacion de integridad de la poblacion
        if len(prev_step_data['pos']) != len(current_step_data['pos']):
            print(f"ADVERTENCIA: Cambio en el numero de particulas ({len(prev_step_data['pos'])} -> {len(current_step_data['pos'])})")
            # Podria ocurrir si hay salida de particulas del dominio
            prev_step_data = current_step_data
            continue

        print(f"Calculando metricas de movimiento: {os.path.basename(step_files[i-1])} -> {os.path.basename(current_step_path)}")

        # Calcular desplazamientos (Euclideo)
        distances = np.linalg.norm(current_step_data['pos'] - prev_step_data['pos'], axis=1)

        mean_disp = np.mean(distances)
        max_disp = np.max(distances)
        anomalies = np.sum(distances > MAX_ALLOWED_DISPLACEMENT)

        print(f"  RESULTADOS:")
        print(f"    Desplazamiento medio: {mean_disp:.6e}")
        print(f"    Desplazamiento max:   {max_disp:.6e}")

        if anomalies > 0:
            print(f"    ALERTA: {anomalies} particulas superan el umbral de movimiento de {MAX_ALLOWED_DISPLACEMENT}")
        else:
            print(f"    Estado: Movimiento fluido y consistente.")

        print("-" * 50)

        # El actual pasa a ser el previo para la siguiente iteracion
        prev_step_data = current_step_data

if __name__ == "__main__":
    start_time = time.time()
    analyze_motion()
    end_time = time.time()
    print(f"PROCESO FINALIZADO")
    print(f"Tiempo total de ejecucion: {(end_time - start_time)/60:.2f} minutos")