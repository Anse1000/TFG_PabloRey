# TFG_PRA — Trabajo de Fin de Grado Pablo Rey Ansemil

## Estructura del repositorio

- `src/` — Código fuente principal (C/CUDA)
     - `src/main.c` — Punto de entrada del programa
     - `src/cpu/` — Implementación CPU (p. ej., estructuras tipo *octree* y simulación)
     - `src/gpu/` — Implementación GPU (CUDA: `.cu`, `.cuh`)
     - `src/*.c, src/*.h` — Utilidades comunes: lectura/escritura, cálculos, tipos, etc.
- `scripts/` — Scripts auxiliares (Python y Bash) para preparar datos, comparar resultados y automatizar pruebas/render
- `docs/` — Memoria y materiales en LaTeX (incluye `memoria_tfg.tex` y subcarpetas de contenido)
- `pruebas/` / `test_results/` — Pruebas y salidas de validación/benchmark
- `CMakeLists.txt` — Configuración de compilación del proyecto
- `*.slurm` — Lanzadores para ejecución en clúster (Slurm)

---

## Requisitos
- **CMake** (recomendado >= 3.20)
- **Librería HDF5** (headers y librería enlazable)
  - Si la instalas en una ruta no estándar, configura `HDF5_ROOT` o `CMAKE_PREFIX_PATH`.
- **OpenMP** (soporte del compilador; en GCC suele venir integrado, en Intel también)
### Para compilar (CPU)
- **Compilador Intel (icc o icx)** recomendado
> **Nota:** si `icc` no está disponible, el proyecto puede compilarse con **GCC (`gcc`)** (manteniendo OpenMP y HDF5 configurados).
### Para compilar con GPU (CUDA)
- **NVIDIA CUDA Toolkit** (con `nvcc`)
- Drivers NVIDIA actualizados y GPU compatible
### Utilidades opcionales
- **Python 3** para ejecutar herramientas en `scripts/`
---

## Compilación

El proyecto usa la opción de CMake `ENABLE_CUDA` para elegir entre implementaciones CPU y GPU:

- **CPU**: compilar con `-DENABLE_CUDA=OFF`
- **GPU (CUDA)**: compilar con `-DENABLE_CUDA=ON`

Ejemplo:
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=OFF -S . -B ./build
cmake --build ./build
``` 
## Ejecución
Uso:
```bash
./build/TFG_PRA <archivo_estrellas> <archivo_salida> <STEPS>
``` 
Parámetros:
- `<archivo_estrellas>`: **carpeta** con los datos de entrada de las estrellas.
- `<archivo_salida>`: **carpeta** donde se guardarán los resultados.
- `<STEPS>`: número de pasos de la simulación (entero > 0).

## Documentación (memoria)

La memoria del TFG se encuentra en `docs/` (LaTeX). Se compila desde el archivo principal:
- `docs/memoria_tfg.tex`

---