#include "cuda_functions.cuh"
#include "aux_fun.h"
#include <sys/time.h>

double MIN_NODE_SIZE_UNIFIED = 1e-10;

#define BLOCK_SIZE 256

long octree_new_node_unified(Octree *tree, float cx, float cy, float cz, float half_size) {
    if (tree->size >= tree->capacity) {
        fprintf(stderr, "Error: No hay suficiente memoria para crear un nuevo nodo\n");
        exit(1);
    }
    size_t i = tree->size++;

    tree->center_x[i] = cx;
    tree->center_y[i] = cy;
    tree->center_z[i] = cz;
    tree->half_size[i] = half_size;
    tree->mass[i] = 0.0F;
    tree->com_x[i] = 0.0;
    tree->com_y[i] = 0.0;
    tree->com_z[i] = 0.0;
    tree->star_index[i] = -1;

    for (int j = 0; j < 8; j++)
        tree->children[i][j] = INVALID_INDEX;
    return i;
}

void octree_insert_unified(Octree *tree, Star *stars, long node_index, long star_index) {
    float cx = tree->center_x[node_index];
    float cy = tree->center_y[node_index];
    float cz = tree->center_z[node_index];
    float hs = tree->half_size[node_index];

    double x = stars->Cx[star_index];
    double y = stars->Cy[star_index];
    double z = stars->Cz[star_index];
    float m = stars->mass[star_index];

    float old_mass = tree->mass[node_index];
    float new_mass = old_mass + m;

    tree->com_x[node_index] = (tree->com_x[node_index] * old_mass + x * m) / new_mass;
    tree->com_y[node_index] = (tree->com_y[node_index] * old_mass + y * m) / new_mass;
    tree->com_z[node_index] = (tree->com_z[node_index] * old_mass + z * m) / new_mass;
    tree->mass[node_index] = new_mass;

    if (hs * 2.0 <= MIN_NODE_SIZE_UNIFIED) {
        if (tree->star_index[node_index] == -1)
            tree->star_index[node_index] = star_index;
        return;
    }

    int oct = ((x >= cx) << 2) | ((y >= cy) << 1) | (z >= cz);

    if (tree->children[node_index][oct] == INVALID_INDEX) {
        float offset = hs * 0.5F;
        float new_cx = cx + ((oct & 4) ? offset : -offset);
        float new_cy = cy + ((oct & 2) ? offset : -offset);
        float new_cz = cz + ((oct & 1) ? offset : -offset);

        long child_index = octree_new_node_unified(tree, new_cx, new_cy, new_cz, offset);
        tree->children[node_index][oct] = child_index;

        tree->star_index[child_index] = star_index;
        tree->mass[child_index] = m;
        tree->com_x[child_index] = x;
        tree->com_y[child_index] = y;
        tree->com_z[child_index] = z;
    } else {
        long child = tree->children[node_index][oct];
        if (tree->star_index[child] >= 0) {
            long existing_star = tree->star_index[child];
            tree->star_index[child] = -1;

            tree->mass[child] = 0.0F;
            tree->com_x[child] = 0.0;
            tree->com_y[child] = 0.0;
            tree->com_z[child] = 0.0;

            octree_insert_unified(tree, stars, child, existing_star);
            octree_insert_unified(tree, stars, child, star_index);
        } else {
            octree_insert_unified(tree, stars, child, star_index);
        }
    }
}
// Construir árbol optimizado para GPU
Octree *build_tree_unified_memory(Star *stars) {
    struct timeval start, end;

    gettimeofday(&start, NULL);

    int device_count;
    cudaGetDeviceCount(&device_count);

    printf("Construyendo árbol en Unified memory\n");
    fflush(stdout);

    // Crear árbol en unified memory
    Octree *tree;
    cudaMallocManaged(&tree, sizeof(Octree));
    memset(tree, 0, sizeof(Octree));

    //Asignar una capacidad inicial que tenga sentido
    tree->capacity = stars->size*1.5f;
    tree->size = 0;

    cudaMallocManaged(&tree->center_x, tree->capacity * sizeof(float));
    cudaMallocManaged(&tree->center_y, tree->capacity * sizeof(float));
    cudaMallocManaged(&tree->center_z, tree->capacity * sizeof(float));
    cudaMallocManaged(&tree->half_size, tree->capacity * sizeof(float));
    cudaMallocManaged(&tree->mass, tree->capacity * sizeof(float));
    cudaMallocManaged(&tree->com_x, tree->capacity * sizeof(double));
    cudaMallocManaged(&tree->com_y, tree->capacity * sizeof(double));
    cudaMallocManaged(&tree->com_z, tree->capacity * sizeof(double));
    cudaMallocManaged(&tree->children, tree->capacity * sizeof(unsigned int[8]));
    cudaMallocManaged(&tree->star_index, tree->capacity * sizeof(long));

    // Inicializar
    for (size_t i = 0; i < tree->capacity; i++) {
        for (int j = 0; j < 8; j++) tree->children[i][j] = INVALID_INDEX;
        tree->star_index[i] = -1;
    }

    // Construir árbol usando funciones existentes
    float cx, cy, cz, hs;
    compute_root_bounds(stars, &cx, &cy, &cz, &hs, &MIN_NODE_SIZE_UNIFIED);
    printf("MIN_NODE_SIZE_UNIFIED: %f\n", MIN_NODE_SIZE_UNIFIED); fflush(stdout);

    long root = octree_new_node_unified(tree, cx, cy, cz, hs);

    for (unsigned long i = 0; i < stars->size; i++) {
        octree_insert_unified(tree, stars, root, i);
        if (i%1000000==0) printf("Insertado %ld estrellas\n",i); fflush(stdout);
    }

    // Configurar hints para acceso desde todas las GPUs
    for (int i = 0; i < device_count; i++) {
        // Prefetch y configurar políticas de acceso
        cudaMemPrefetchAsync(tree, sizeof(Octree), i);
        cudaMemPrefetchAsync(tree->center_x, tree->capacity * sizeof(float), i);
        cudaMemPrefetchAsync(tree->center_y, tree->capacity * sizeof(float), i);
        cudaMemPrefetchAsync(tree->center_z, tree->capacity * sizeof(float), i);
        cudaMemPrefetchAsync(tree->half_size, tree->capacity * sizeof(float), i);
        cudaMemPrefetchAsync(tree->mass, tree->capacity * sizeof(float), i);
        cudaMemPrefetchAsync(tree->com_x, tree->capacity * sizeof(double), i);
        cudaMemPrefetchAsync(tree->com_y, tree->capacity * sizeof(double), i);
        cudaMemPrefetchAsync(tree->com_z, tree->capacity * sizeof(double), i);
        cudaMemPrefetchAsync(tree->children, tree->capacity * sizeof(unsigned int[8]), i);
        cudaMemPrefetchAsync(tree->star_index, tree->capacity * sizeof(long), i);

        // Configurar como solo lectura para optimizar caching
        cudaMemAdvise(tree->center_x, tree->capacity * sizeof(float), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->center_y, tree->capacity * sizeof(float), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->center_z, tree->capacity * sizeof(float), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->half_size, tree->capacity * sizeof(float), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->mass, tree->capacity * sizeof(float), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->com_x, tree->capacity * sizeof(double), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->com_y, tree->capacity * sizeof(double), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->com_z, tree->capacity * sizeof(double), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->children, tree->capacity * sizeof(unsigned int[8]), cudaMemAdviseSetReadMostly, i);
        cudaMemAdvise(tree->star_index, tree->capacity * sizeof(long), cudaMemAdviseSetReadMostly, i);
    }

    gettimeofday(&end, NULL);
    double secs = get_seconds(start, end);
    size_t actual_memory = tree->size * (sizeof(float) * 4 + sizeof(float) + sizeof(double) * 3 + sizeof(unsigned int) *
                                         8 + sizeof(long));

    printf("Árbol construido: %ld nodos en %.2f segundos (%.1f GB)\n",
           tree->size, secs, actual_memory / (1024.0 * 1024.0 * 1024.0));
    fflush(stdout);

    return tree;
}

__global__ void compute_forces_kernel(const Octree *tree, size_t star_count, 
                                      double *Cx, double *Cy, double *Cz,
                                      double *ax, double *ay, double *az,
                                      double theta) {
    int star_idx = blockIdx.x * blockDim.x + threadIdx.x;
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

void free_tree_unified(Octree *tree) {
    if (tree != NULL) {
        if (tree->center_x != NULL) {
            cudaFree(tree->center_x);
            cudaFree(tree->center_y);
            cudaFree(tree->center_z);
            cudaFree(tree->half_size);
            cudaFree(tree->mass);
            cudaFree(tree->com_x);
            cudaFree(tree->com_y);
            cudaFree(tree->com_z);
            cudaFree(tree->children);
            cudaFree(tree->star_index);
        }
        cudaFree(tree);
    }
}

//reordenar estrellas por octante para hacer calculos en gpu
void reorder_stars(Star *stars, Octree *tree, unsigned int *offsets) {
    float cx = tree->center_x[0];
    float cy = tree->center_y[0];
    float cz = tree->center_z[0];

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

__host__ void compute_acceleration_multi_gpu(long N, int iterations, double *ax, double *ay, double *az,
                                             const unsigned int *offsets, int device_count, Octree *tree,
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
            int grid_size = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
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
extern "C" void simulate_multi_gpu_unified(Star *estrellas, const long N, const char *outputfile) {
    struct timeval start, end;

    gettimeofday(&start, NULL);
    int device_count = 0;
    cudaGetDeviceCount(&device_count);

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

    for (int step = 0; step < STEPS; step++) {
        printf("\n--- Paso %d ---\n", step + 1);
        // Construir árbol
        Octree *octree = build_tree_unified_memory(estrellas);
        unsigned int offsets[8];
        reorder_stars(estrellas, octree, offsets);
        printf("Iniciando fase 1\n"); fflush(stdout);
        compute_acceleration_multi_gpu(N, iterations, ax, ay, az, offsets, device_count, octree, estrellas, streams);
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
        free_tree_unified(octree);
        octree = build_tree_unified_memory(estrellas);
        printf("Iniciando fase 2\n"); fflush(stdout);
        compute_acceleration_multi_gpu(N, iterations, ax, ay, az, offsets, device_count, octree, estrellas, streams);
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
    for (int i = 0; i < N; i++) {
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n", estrellas->id[i], estrellas->Cx[i], estrellas->Cy[i],
                estrellas->Cz[i]);
    }
}