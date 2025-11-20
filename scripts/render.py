import os
import glob
from paraview.simple import *
import argparse

GAUSSIAN_RADIUS = 0.0010   # radio base de Gaussian Points
COLOR_FIELD = "MASS"        # columna para color y escala
RES_X = 1920
RES_Y = 1080

def process_step(step_dir):
    # Cargar todos los CSV del step
    csv_files = sorted(glob.glob(os.path.join(step_dir, "*.csv")))
    if not csv_files:
        return None

    point_sources = []
    for csv_file in csv_files:
        r = CSVReader(FileName=csv_file)
        tableToPoints = TableToPoints(Input=r)
        tableToPoints.XColumn = "X"
        tableToPoints.YColumn = "Y"
        tableToPoints.ZColumn = "Z"
        point_sources.append(tableToPoints)

    if len(point_sources) == 1:
        merged = point_sources[0]
    else:
        merged = AppendDatasets(Input=point_sources)

    return merged

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", required=True, help="Carpeta que contiene step_0, step_1, ...")
    parser.add_argument("--outdir", required=True, help="Carpeta donde guardar imágenes PNG")
    parser.add_argument("--resx", type=int, default=RES_X)
    parser.add_argument("--resy", type=int, default=RES_Y)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Detectar steps
    step_dirs = sorted([os.path.join(args.indir, d) for d in os.listdir(args.indir) if d.startswith("step_")])
    if not step_dirs:
        print("No se encontraron steps")
        exit(1)

    # Crear vista
    view = CreateView("RenderView")
    view.ViewSize = [args.resx, args.resy]
    view.UseColorPaletteForBackground = 0
    view.Background = [0.0, 0.0, 0.0]

    # --- Calcular bounding box global ---
    global_bounds = [1e99, -1e99, 1e99, -1e99, 1e99, -1e99]
    for step_dir in step_dirs:
        dataset = process_step(step_dir)
        if dataset is None:
            continue
        info = dataset.GetDataInformation().GetBounds()
        for j in range(3):
            global_bounds[2*j]   = min(global_bounds[2*j], info[2*j])
            global_bounds[2*j+1] = max(global_bounds[2*j+1], info[2*j+1])
        Delete(dataset)
        del dataset

    # --- Configurar cámara al bounding global ---
    cx = 0.5 * (global_bounds[0] + global_bounds[1])
    cy = 0.5 * (global_bounds[2] + global_bounds[3])
    cz = 0.5 * (global_bounds[4] + global_bounds[5])
    dx = global_bounds[1] - global_bounds[0]
    dy = global_bounds[3] - global_bounds[2]
    dz = global_bounds[5] - global_bounds[4]
    dist = max(dx, dy, dz)

    view.CameraFocalPoint = [cx, cy, cz]
    view.CameraPosition   = [cx, cy, cz + dist]
    view.CameraViewUp     = [0, 1, 0]

    # --- Render paso a paso ---
    for i, step_dir in enumerate(step_dirs):
        print(f"Procesando {step_dir} ...")
        dataset = process_step(step_dir)
        if dataset is None:
            continue

        display = Show(dataset, view, 'GeometryRepresentation')
        display.Representation = 'Point Gaussian'
        display.GaussianRadius = GAUSSIAN_RADIUS
        display.Emissive = 1
        display.ScaleByArray = 1
        display.SetScaleArray = ['POINTS', COLOR_FIELD]
        display.UseScaleFunction = 0

        ColorBy(display, ('POINTS', COLOR_FIELD))
        lut = GetColorTransferFunction(COLOR_FIELD)
        lut.ApplyPreset("Black-Body Radiation", True)

        # --- Escala logarítmica ---
        lut.MapControlPointsToLogSpace()
        lut.UseLogScale = 1

        # --- Ajustar rango para evitar "todo rojo" ---
        lut.RescaleTransferFunction(0.05, 5)

        display.LookupTable = lut
        display.SetScalarBarVisibility(view, True)

        Render()
        filename = os.path.join(args.outdir, f"step_{i:04d}.png")
        SaveScreenshot(filename, view, ImageResolution=[args.resx, args.resy])
        print(f" -> Guardada {filename}")

        # --- Liberar memoria ---
        Hide(dataset, view)
        Delete(dataset)
        del dataset

    print("Renderizado por pasos completado")
