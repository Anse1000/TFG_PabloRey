import sys
import pandas as pd
import numpy as np

# Archivos de entrada
file1, file2, outfile = sys.argv[1:]

# Leer CSV
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Asegurar orden por ID
df1 = df1.sort_values("ID").reset_index(drop=True)
df2 = df2.sort_values("ID").reset_index(drop=True)

# Verificación de IDs
if not (df1["ID"].equals(df2["ID"])):
    raise ValueError("IDs no coinciden entre los dos CSV.")

# Calcular diferencias
diffs = df1[["X", "Y", "Z"]].values - df2[["X", "Y", "Z"]].values
df_diff = pd.DataFrame(diffs, columns=["dx", "dy", "dz"])

# Distancia euclídea entre coordenadas
df_diff["dist"] = np.linalg.norm(diffs, axis=1)

# Función auxiliar para calcular métricas con percentiles
def metrics(series, name):
    percentiles = [1, 5, 25, 50, 75, 95, 99]
    stats = {
        f"{name}_mean": np.mean(series),
        f"{name}_std": np.std(series),
        f"{name}_min": np.min(series),
        f"{name}_max": np.max(series),
        f"{name}_rmse": np.sqrt(np.mean(series**2)),
    }
    for p in percentiles:
        stats[f"{name}_p{p}"] = np.percentile(series, p)
    return stats

# Calcular métricas para cada dimensión y la distancia total
all_stats = {}
for col in ["dx", "dy", "dz", "dist"]:
    all_stats.update(metrics(df_diff[col], col))

# Guardar en CSV
pd.DataFrame([all_stats]).to_csv(outfile, index=False)
