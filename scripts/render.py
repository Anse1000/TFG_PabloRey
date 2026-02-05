import os
import sys
import argparse
import gc
from paraview.simple import *
from paraview.modules.vtkRemotingCore import vtkProcessModule

# --- UTILIDADES MPI ---
pm = vtkProcessModule.GetProcessModule()
rank = pm.GetPartitionId()

def print0(msg):
    if rank == 0:
        print(f"[Rank 0] {msg}", flush=True)

# --- CONFIGURACIÓN ---
GAUSSIAN_BASE_RADIUS = 0.005
COLOR_FIELD = "MASS"
RES_X = 15360 # Usando tus resoluciones altas
RES_Y = 8640

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--resx", type=int, default=RES_X)
    parser.add_argument("--resy", type=int, default=RES_Y)

    try:
        args = parser.parse_args()
    except SystemExit:
        if rank == 0: raise
        sys.exit(0)

    # 1. Guardamos rutas absolutas (CRÍTICO)
    abs_indir = os.path.abspath(args.indir)
    abs_outdir = os.path.abspath(args.outdir)
    original_cwd = os.getcwd() # Guardamos donde empezamos

    if rank == 0 and not os.path.exists(abs_outdir):
        os.makedirs(abs_outdir)

    # Listamos carpetas desde la ruta absoluta
    all_files = sorted([d for d in os.listdir(abs_indir) if d.startswith("step_")])

    print0(f"Procesando {len(all_files)} pasos (Estrategia: Cambiar Directorio)...")

    # Crear vista una sola vez
    view = CreateView("RenderView")
    view.ViewSize = [args.resx, args.resy]
    view.Background = [0.0, 0.0, 0.0]
    view.OrientationAxesVisibility = 0
    view.RemoteRenderThreshold = 0 # Fuerza a que no se envíen datos entre nodos
    view.UseCache = 0      # Desactiva copias de geometría de baja resolución   

    for step_dir_name in all_files:
        try:
            step_id = int(step_dir_name.split('_')[-1])
        except ValueError: continue

        # Definimos rutas
        step_folder_path = os.path.join(abs_indir, step_dir_name)
        xmf_filename = f"step_{step_id:04d}.xmf"
        png_file = os.path.join(abs_outdir, f"step_{step_id:04d}.png")

        # Verificar existencia
        full_xmf_path = os.path.join(step_folder_path, xmf_filename)
        if not os.path.exists(full_xmf_path): continue
        if os.path.exists(png_file):
            print0(f"Saltando {step_id}, ya existe.")
            continue

        try:
            # --- TRUCO CRÍTICO: CAMBIAR EL DIRECTORIO DE TRABAJO ---
            # Nos movemos a la carpeta donde están los datos.
            # Así el lector encuentra los .h5 locales sin problemas de rutas.
            os.chdir(step_folder_path)

            # --- MODIFICACIÓN DEL READER ---
            # Intentamos usar Xdmf3ReaderS (Correcto para ParaView 5.13+ y XDMF 3.0)
            # Si falla (versiones viejas), usamos el Legacy.
            try:
                # Nota: Xdmf3ReaderS usa 'FileName' (singular) y string directo
                reader = Xdmf3ReaderS(FileName=full_xmf_path)
            except NameError:
                # Fallback: XDMFReader usa 'FileNames' (plural) y lista
                reader = XDMFReader(FileNames=[full_xmf_path])

            reader.UpdatePipeline()

            # --- RENDERIZADO ---
            display = Show(reader, view, 'GeometryRepresentation')
            display.Representation = 'Point Gaussian'
            display.SetScaleArray = ['POINTS', COLOR_FIELD]
            display.GaussianRadius = GAUSSIAN_BASE_RADIUS

            ColorBy(display, ('POINTS', COLOR_FIELD))
            lut = GetColorTransferFunction(COLOR_FIELD)
            lut.ApplyPreset('Black-Body Radiation', 1)
            lut.RescaleTransferFunction(0.3, 6.0)
            display.SetScalarBarVisibility(view, False)

            # Cámara (Ajustar según tus necesidades)
            view.CameraFocalPoint = [0, 0, 0]
            view.ResetCamera()
            # Forzamos posición 'lejana' si ResetCamera falla o queda muy cerca
            # view.CameraPosition = [0, 0, 10.0]

            Render()

            if rank == 0:
                SaveScreenshot(png_file, view, ImageResolution=[args.resx, args.resy])
                print(f" -> Guardado: {png_file}")

            Delete(display)
            Delete(reader)

            # Volvemos al directorio original por seguridad antes del siguiente paso
            os.chdir(original_cwd)

            # Limpieza
            gc.collect()

        except Exception as e:
            # Aseguramos volver al directorio original si falla
            os.chdir(original_cwd)
            sys.stderr.write(f"ERROR {step_id}: {e}\n")

    print0("Fin del proceso.")