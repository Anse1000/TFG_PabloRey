#include "cuda_functions.cuh"
#include "aux_fun.h"
#include <sys/time.h>
#include "octree.h"

#define BLOCK_SIZE 256

__global__ void compute_forces_kernel(const Octree *tree, size_t star_count, 
                                      const double *Cx, const double *Cy, const double *Cz,
                                      double *ax, double *ay, double *az,
                                      double theta) {
    unsigned int star_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (star_idx >= star_count) return;

    // Inicializar aceleraciones
    double acc_x = 0.0, acc_y = 0.0, acc_z = 0.0;
    
    // Stack explícito para la traversal del árbol
    const int MAX_STACK_SIZE = 64;
    long stack[MAX_STACK_SIZE];
    int top = -1;
    
    // Inicializar con la raíz
    stack[++top] = 0;

    while (top >= 0) {
        long node_idx = stack[top--];

        // Saltar si es la misma estrella
        if (tree->star_index[node_idx] == star_idx)
            continue;

        double dx = tree->com_x[node_idx] - Cx[star_idx];
        double dy = tree->com_y[node_idx] - Cy[star_idx];
        double dz = tree->com_z[node_idx] - Cz[star_idx];
        double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;
        double dist = sqrt(dist_sq);

        double s = 2.0 * tree->half_size[node_idx];

        // Aplicar criterio de Barnes-Hut
        if ((s / dist) < theta || tree->star_index[node_idx] >= 0) {
            // Calcular fuerza
            double force = -G * tree->mass[node_idx] / (dist_sq * dist);
            acc_x += force * dx;
            acc_y += force * dy;
            acc_z += force * dz;
        } else {
            // Añadir hijos al stack
            for (int i = 0; i < 8; i++) {
                unsigned int child = tree->children[node_idx][i];
                if (child != INVALID_INDEX && top + 1 < MAX_STACK_SIZE) {
                    stack[++top] = child;
                }
            }
        }
    }
    
    // Escribir resultados
    ax[star_idx] = acc_x;
    ay[star_idx] = acc_y;
    az[star_idx] = acc_z;
}

inline int get_octant(double cx, double cy, double cz, double x, double y, double z) {
    return ((x >= cx) << 2) | ((y >= cy) << 1) | (z >= cz);
}

//reordenar estrellas por octante para hacer calculos en gpu
void reorder_stars(Star *stars,float cx, float cy, float cz, unsigned int *offsets) {
    size_t counts[8] = {0};
    for (size_t i = 0; i < stars->size; i++) {
        int oct = get_octant(cx, cy, cz, stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        counts[oct]++;
    }
    offsets[0] = 0;
    for (int i = 1; i < 8; i++) {
        offsets[i] = offsets[i - 1] + counts[i - 1];
    }
    unsigned int ends[8];
    memcpy(ends, offsets, 8*sizeof(unsigned int));
    for (size_t i = 0; i < stars->size;) {
        int oct = get_octant(cx, cy, cz, stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        if (i >= offsets[oct] && i < ends[oct]) {
            // Ya está en su rango
            i++;
        } else {
            // Debe ir en ends[oct]
            size_t dest = ends[oct];
            swap_star_elements(stars, i, dest);
            ends[oct]++;
        }
    }
    printf("Estrellas ordenadas por octante\n"); fflush(stdout);
}

__host__ void compute_acceleration_multi_gpu(unsigned int N, int iterations, double *ax, double *ay, double *az,
                                             const unsigned int *offsets, int device_count, Octree **trees,
                                             Star *estrellas, cudaStream_t *streams) {
    for (int i = 0; i < iterations; i++) {
        for (int dev = 0; dev < device_count; dev++) {
            int octant = i * device_count + dev;
            if (octant >= 8) break;
            
            cudaSetDevice(dev);

            long start = offsets[octant];
            long end = (octant == 7) ? N : offsets[octant + 1];
            long count = end - start;

            if (count == 0) continue;

            // Reservar memoria en GPU
            Octree *tree = trees[octant];

            double *d_x, *d_y, *d_z;
            double *d_ax, *d_ay, *d_az;
            cudaMalloc(&d_x, count * sizeof(double));
            cudaMalloc(&d_y, count * sizeof(double));
            cudaMalloc(&d_z, count * sizeof(double));
            cudaMalloc(&d_ax, count * sizeof(double));
            cudaMalloc(&d_ay, count * sizeof(double));
            cudaMalloc(&d_az, count * sizeof(double));


            printf("Memoria reservada en GPU it: %d dev: %d count: %ld\n",i,dev,count); fflush(stdout);

            // Copiar datos a GPU
            cudaMemcpyAsync(d_x, &estrellas->Cx[start], count * sizeof(double), 
                           cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_y, &estrellas->Cy[start], count * sizeof(double), 
                           cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_z, &estrellas->Cz[start], count * sizeof(double), 
                           cudaMemcpyHostToDevice, streams[dev]);

            // Inicializar aceleraciones a cero
            cudaMemsetAsync(d_ax, 0, count * sizeof(double), streams[dev]);
            cudaMemsetAsync(d_ay, 0, count * sizeof(double), streams[dev]);
            cudaMemsetAsync(d_az, 0, count * sizeof(double), streams[dev]);

            // Esperar a que termine la copia antes de lanzar kernel
            cudaStreamWaitEvent(streams[dev], 0, 0);
            printf("Copia terminada en GPU it:%d dev:%d\n",i,dev); fflush(stdout);

            // Lanzar kernel
            unsigned int grid_size = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
            compute_forces_kernel<<<grid_size, BLOCK_SIZE, 0, streams[dev]>>>(
                tree, count, d_x, d_y, d_z, d_ax, d_ay, d_az, 0.2);

            // Verificar errores del kernel
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("Error en kernel: %s\n", cudaGetErrorString(err));
            }

            // Copiar resultados
            cudaMemcpyAsync(&ax[start], d_ax, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);
            cudaMemcpyAsync(&ay[start], d_ay, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);
            cudaMemcpyAsync(&az[start], d_az, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);

            cudaStreamSynchronize(streams[dev]);

            printf("Copia resultados terminada it:%d dev:%d\n",i,dev); fflush(stdout);

            // Liberar memoria
            cudaFree(d_x);
            cudaFree(d_y);
            cudaFree(d_z);
            cudaFree(d_ax);
            cudaFree(d_ay);
            cudaFree(d_az);
        }
        
        // Sincronizar todos los dispositivos
        for (int dev = 0; dev < device_count; dev++) {
            cudaSetDevice(dev);
            cudaStreamSynchronize(streams[dev]);
        }
    }
}

// Función principal de simulacion en gpus
extern "C" void simulate_multi_gpu_unified(Star *estrellas,const int steps, const long N, const char *outputfile) {
    struct timeval start, end;
    float cx,cy,cz;
    gettimeofday(&start, NULL);
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    if (device_count==0) {
        printf("No hay dispositivos disponibles\n");
        exit(1);
    }

    // Arrays para resultados
    auto *ax = static_cast<double *>(malloc(N * sizeof(double)));
    auto *ay = static_cast<double *>(malloc(N * sizeof(double)));
    auto *az = static_cast<double *>(malloc(N * sizeof(double)));

    // Crear streams para cada GPU
    auto *streams = static_cast<cudaStream_t *>(malloc(device_count * sizeof(cudaStream_t)));
    for (int i = 0; i < device_count; i++) {
        cudaSetDevice(i);
        cudaStreamCreate(&streams[i]);
    }
    const int iterations = (8 + device_count - 1) / device_count;

    printf("=== Iniciando simulación con %d GPUs ===\n", device_count);

    for (int step = 0; step < steps; step++) {
        printf("\n--- Paso %d ---\n", step + 1);
        // Construir árbol
        Octree **octrees = build_tree_gpu(estrellas,&cx,&cy,&cz);
        unsigned int offsets[8];
        reorder_stars(estrellas, cx,cy,cz, offsets);
        printf("Iniciando fase 1\n"); fflush(stdout);
        compute_acceleration_multi_gpu(N, iterations, ax, ay, az, offsets, device_count, octrees, estrellas, streams);
        // Aplicar integración Leapfrog
        double DT2 = 0.5 * DT;
        for (long i = 0; i < N; i++) {
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);

            estrellas->Cx[i] = fma(DT, estrellas->Vx[i], estrellas->Cx[i]);
            estrellas->Cy[i] = fma(DT, estrellas->Vy[i], estrellas->Cy[i]);
            estrellas->Cz[i] = fma(DT, estrellas->Vz[i], estrellas->Cz[i]);
        }
        octrees = build_tree_gpu(estrellas,&cx,&cy,&cz);
        reorder_stars(estrellas, cx,cy,cz, offsets);
        printf("Iniciando fase 2\n"); fflush(stdout);
        compute_acceleration_multi_gpu(N, iterations, ax, ay, az, offsets, device_count, octrees, estrellas, streams);
        for (long i = 0; i < N; i++) {
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
        }
        printf("Paso %d completado\n", step + 1);
        fflush(stdout);
    }
    // Limpiar streams
    for (int dev = 0; dev < device_count; dev++) {
        cudaSetDevice(dev);
        cudaStreamDestroy(streams[dev]);
    }
    gettimeofday(&end, NULL);
    double seconds = get_seconds(start, end);
    int hours = static_cast<int>(seconds / 3600);
    int minutes = (static_cast<int>(seconds) % 3600) / 60;
    double remaining_seconds = fmod(seconds, 60.0);
    printf("Simuladas %ld estrellas en %02d:%02d:%05.2f (hh:mm:ss) usando %d GPUs\n",
           N, hours, minutes, remaining_seconds, device_count);
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
    FILE *file = fopen(outputfile, "w");
    if (!file) {
        fprintf(stderr, "Error al abrir el archivo de salida\n");
        exit(1);
    }
    printf("Escribiendo resultados en %s\n", outputfile); fflush(stdout);
    for (int i = 0; i < N; i++) {
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n", estrellas->id[i], estrellas->Cx[i], estrellas->Cy[i],
                estrellas->Cz[i]);
    }
}