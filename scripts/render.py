import os
import sys
import argparse
from paraview.simple import *
from paraview.modules.vtkRemotingCore import vtkProcessModule

# --- UTILIDADES MPI ---
# Obtener información del proceso MPI actual
pm = vtkProcessModule.GetProcessModule()
rank = pm.GetPartitionId()
nranks = pm.GetNumberOfLocalPartitions()

def print0(msg):
    """Solo imprime si somos el proceso maestro (Rank 0)"""
    if rank == 0:
        print(f"[Rank 0] {msg}")
        sys.stdout.flush()

# --- CONFIGURACIÓN ---
GAUSSIAN_RADIUS = 0.0010
COLOR_FIELD = "MASS"
RES_X = 1920
RES_Y = 1080

if __name__ == "__main__":
    # Argumentos (Parsear solo en Rank 0 y difundir sería ideal,
    # pero simple args funciona si todos reciben los mismos flags)
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", required=True, help="Carpeta con subcarpetas step_XXXX")
    parser.add_argument("--outdir", required=True, help="Carpeta de salida")
    parser.add_argument("--resx", type=int, default=RES_X)
    parser.add_argument("--resy", type=int, default=RES_Y)

    # Para evitar conflictos en MPI al parsear ayuda, envolvemos en try
    try:
        args = parser.parse_args()
    except SystemExit:
        # Evitar que ranks > 0 maten el proceso ruidosamente si falta un arg
        if rank == 0: raise
        sys.exit(0)

    # Rutas absolutas
    args.indir = os.path.abspath(args.indir)
    args.outdir = os.path.abspath(args.outdir)

    if rank == 0:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

    # Listar directorios (asumimos que la estructura de disco es visible para todos)
    # Es seguro hacerlo en todos los ranks si el filesystem es compartido.
    all_files = sorted([d for d in os.listdir(args.indir) if d.startswith("step_")])

    if not all_files:
        print0("No se encontraron carpetas step_XXXX.")
        sys.exit(1)

    # --- SETUP INICIAL DE PARAVIEW ---
    # Crear la vista una sola vez
    view = CreateView("RenderView")
    view.ViewSize = [args.resx, args.resy]
    view.Background = [0.0, 0.0, 0.0]
    view.UseColorPaletteForBackground = 0

    # Optimizaciones para renderizado off-screen
    view.OrientationAxesVisibility = 0

    camera_set = False

    print0(f"Iniciando renderizado de {len(all_files)} pasos con {nranks} procesos MPI...")

    # --- BUCLE DE RENDER ---
    for step_dir_name in all_files:
        step_path = os.path.join(args.indir, step_dir_name)

        # Extraer ID
        try:
            step_id = int(step_dir_name.split('_')[-1])
        except ValueError:
            continue

        xmf_filename = f"step_{step_id:04d}.xmf"
        xmf_filepath = os.path.join(step_path, xmf_filename)
        png_file = os.path.join(args.outdir, f"step_{step_id:04d}.png")

        if not os.path.exists(xmf_filepath):
            continue

        # 1. Cargar Datos
        # ParaView maneja la distribución de datos automáticamente si el formato lo soporta
        reader = XDMFReader(FileNames=[xmf_filepath])

        # 2. Representación Visual
        display = Show(reader, view, 'GeometryRepresentation')
        display.Representation = 'Point Gaussian'
        display.GaussianRadius = GAUSSIAN_RADIUS
        display.Emissive = 1
        display.ScaleByArray = 1
        display.SetScaleArray = ['POINTS', COLOR_FIELD]
        display.UseScaleFunction = 0

        # Color
        ColorBy(display, ('POINTS', COLOR_FIELD))
        lut = GetColorTransferFunction(COLOR_FIELD)
        lut.ApplyPreset("Black-Body Radiation", True)

        # Ajuste Logarítmico
        lut.MapControlPointsToLogSpace()
        lut.UseLogScale = 1
        # Asegúrate de que este rango (0.05, 5.0) tenga sentido para tus datos
        lut.RescaleTransferFunction(0.05, 5.0)

        display.SetScalarBarVisibility(view, True)

        # --- 3. CÁMARA TIPO "VISTA DE GALAXIA" OPTIMIZADA ---

        # El renderizado (Render) final se hará al final del bucle.

        if not camera_set:
            # 1. CRÍTICO: Forzar a ParaView a leer los datos y calcular los bounds
            # Esto es lo que permite a ResetCamera saber dónde están las estrellas.
            reader.UpdatePipeline()

            # 2. Definir A DÓNDE miramos (El centro de tu galaxia)
            view.CameraFocalPoint = [0.0, 0.0, 0.0]

            # 3. Definir EL ÁNGULO (Usamos la vista diagonal elevada)
            view.CameraPosition = [1.0, 1.0, 0.8]
            view.CameraViewUp = [0.0, 0.0, 1.0]

            # 4. Ajustar la distancia automáticamente (Zoom Out)
            # ResetCamera conserva tu ángulo pero ajusta la distancia para que todo quepa.
            view.ResetCamera()

            # 5. (Opcional) Alejar un 20% extra para dar "aire"
            view.GetActiveCamera().Dolly(0.8)

            camera_set = True

            if rank == 0:
                print("  [Cámara] Configuración inicial fijada sin renderizar.")

        # NOTA: Ahora el primer Render() se hará solo después de esta configuración,
        # junto con el resto de la configuración de visualización (colores, LUT, etc.).

        # 4. Render y Guardado
        # Render() sincroniza la imagen entre todos los nodos MPI
        Render()

        if rank == 0:
            SaveScreenshot(png_file, view, ImageResolution=[args.resx, args.resy])
            print(f" -> Guardado: {os.path.basename(png_file)}")

        # 5. Limpieza CRÍTICA
        # Debemos destruir los objetos proxy para liberar memoria antes del siguiente paso
        Delete(display)
        Delete(reader)
        # Opcional: forzar recolección de basura si tienes problemas de RAM
        # import gc; gc.collect()

    print0("Proceso finalizado.")