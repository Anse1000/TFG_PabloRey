from paraview.simple import *
import numpy as np
import sys, os, glob

# --- argumentos desde SLURM ---
if len(sys.argv) < 5:
    print("Uso: pvbatch interpolate_render_mpi.py <dir_inicial> <dir_final> <out_dir> <nframes>")
    sys.exit(1)

dir_init, dir_final, out_dir, nframes = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
os.makedirs(out_dir, exist_ok=True)

# --- buscar chunks ---
chunks_init = sorted(glob.glob(os.path.join(dir_init, "chunk_*.csv")))
chunks_final = sorted(glob.glob(os.path.join(dir_final, "chunk_*.csv")))

if len(chunks_init) != len(chunks_final):
    raise RuntimeError("Número de chunks inicial/final no coincide")

nchunks = len(chunks_init)
rank = int(os.getenv("PMI_RANK", os.getenv("OMPI_COMM_WORLD_RANK", "0")))
size = int(os.getenv("PMI_SIZE", os.getenv("OMPI_COMM_WORLD_SIZE", "1")))

# --- asignar chunks a cada proceso ---
my_chunks = [i for i in range(nchunks) if i % size == rank]

print(f"[RANK {rank}] procesando chunks {my_chunks}")

# --- configuración global de render ---
view = CreateView("RenderView")
view.ViewSize = [1920, 1080]
view.Background = [0, 0, 0]

view.EnableRayTracing = 1
view.BackEnd = 'OSPRay pathtracer'
view.SamplesPerPixel = 8
view.Denoise = 1

camera = GetActiveCamera()
camera.SetPosition(0, -5000, 2000)
camera.SetFocalPoint(0, 0, 0)
camera.SetViewUp(0, 0, 1)

# --- bucle sobre frames ---
for f in range(nframes+1):
    t = f / nframes
    sources = []

    for idx in my_chunks:
        # cargar csv inicial/final
        data0 = np.loadtxt(chunks_init[idx], delimiter=",", skiprows=1)
        data1 = np.loadtxt(chunks_final[idx], delimiter=",", skiprows=1)

        # data formato: ID, X, Y, Z
        interp = data0.copy()
        interp[:,1:4] = (1-t)*data0[:,1:4] + t*data1[:,1:4]

        # guardar temporal para que ParaView lo lea
        tmpfile = f"{out_dir}/rank{rank}_chunk{idx}_frame{f:04d}.csv"
        np.savetxt(tmpfile, interp, delimiter=",", header="ID,X,Y,Z", comments='')

        # crear fuente CSV en ParaView
        src = CSVReader(FileName=[tmpfile])
        tableToPoints = TableToPoints(Input=src)
        tableToPoints.XColumn = "X"
        tableToPoints.YColumn = "Y"
        tableToPoints.ZColumn = "Z"

        glyph = Glyph(Input=tableToPoints, GlyphType="Sphere")
        glyph.GlyphMode = 'All Points'
        glyph.ScaleFactor = 0.5

        Show(glyph, view)
        sources.append((src, tableToPoints, glyph))

    # render y guardar
    view.ResetCamera()
    filename = os.path.join(out_dir, f"frame_{f:04d}.png")
    SaveScreenshot(filename, view, ImageResolution=[1920,1080])

    # limpiar objetos de este frame
    for src, table, glyph in sources:
        Delete(glyph); del glyph
        Delete(table); del table
        Delete(src); del src
