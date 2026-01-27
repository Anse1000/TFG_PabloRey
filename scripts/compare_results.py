import h5py
import numpy as np
import xml.etree.ElementTree as ET
import os
import time
from concurrent.futures import ThreadPoolExecutor

# CONFIGURACIÓN DE RUTAS
PATH_CPU_XMF = "/mnt/lustre/scratch/nlsas/home/ulc/es/pra/resultados/galactic_positions_cpu/step_0001/step_0001.xmf"
PATH_CUDA_XMF = "/mnt/lustre/scratch/nlsas/home/ulc/es/pra/resultados/galactic_positions_cuda/step_0001/step_0001.xmf"
MAX_WORKERS = 64

def get_chunk_info(xdmf_path):
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

def read_h5_chunk(chunk, out_ids, out_pos, out_mass):
    m = chunk['mapping']
    start, end = chunk['start'], chunk['start'] + chunk['n']
    with h5py.File(m['ID'][0], 'r') as f:
        out_ids[start:end] = f[m['ID'][1]][:]
        out_pos[start:end] = f[m['XYZ'][1]][:]
        out_mass[start:end] = f[m['MASS'][1]][:]

def load_dataset_parallel(xdmf_path):
    print(f"--- Cargando: {xdmf_path} ---", flush=True)
    chunks, total_n = get_chunk_info(xdmf_path)
    ids = np.empty(total_n, dtype=np.uint64)
    pos = np.empty((total_n, 3), dtype=np.float64)
    mass = np.empty(total_n, dtype=np.float32)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for chunk in chunks:
            executor.submit(read_h5_chunk, chunk, ids, pos, mass)
    return ids, pos, mass

def run_comparison():
    t_start = time.time()

    # 1. Cargar datos
    id_c, p_c, m_c = load_dataset_parallel(PATH_CPU_XMF)
    id_g, p_g, m_g = load_dataset_parallel(PATH_CUDA_XMF)

    # 2. Alineación RIGUROSA
    print("Alineando datasets por ID...", flush=True)
    idx_c = np.argsort(id_c)
    p_c = p_c[idx_c]; m_c = m_c[idx_c]
    del id_c, idx_c # Liberar

    idx_g = np.argsort(id_g)
    p_g = p_g[idx_g]; m_g = m_g[idx_g]
    del id_g, idx_g

    # 3. Cálculos
    diff = p_g - p_c
    dist_err = np.linalg.norm(diff, axis=1)

    com_c = np.average(p_c, axis=0, weights=m_c)
    com_g = np.average(p_g, axis=0, weights=m_g)
    drift = np.linalg.norm(com_g - com_c)

    r_c = np.linalg.norm(p_c - com_c, axis=1)
    rg_c = np.sqrt(np.sum(m_c * r_c**2) / np.sum(m_c))
    bits = -np.log2(np.median(dist_err) / rg_c)

    print("\n" + "="*60)
    print(f"INFORME (N={len(dist_err):,})")
    print("="*60)
    print(f"MAE Posición:      {np.mean(dist_err):.4e}")
    print(f"Error Máximo:      {np.max(dist_err):.4e}")
    print(f"Drift CoM:    {drift:.4e}")
    print(f"BITS EFECTIVOS:    {bits:.2f}")
    print(f"Tiempo:            {(time.time() - t_start)/60:.2f} min")
    print("="*60)

if __name__ == "__main__":
    run_comparison()