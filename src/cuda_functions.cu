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
    constexpr int MAX_STACK_SIZE = 512;
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

void copy_tree_to_gpu(Octree **d_tree, const Octree *host_tree, cudaStream_t stream) {
    // Verificar que el árbol host no sea nulo
    if (host_tree == NULL) {
        printf("Error: host_tree es NULL\n");
        *d_tree = NULL;
        return;
    }
    
    // 1. Crear estructura en GPU
    Octree *gpu_tree_struct;
    cudaError_t err = cudaMalloc(&gpu_tree_struct, sizeof(Octree));
    if (err != cudaSuccess) {
        printf("Error en cudaMalloc para estructura: %s\n", cudaGetErrorString(err));
        *d_tree = NULL;
        return;
    }

    // 2. Crear una copia temporal en host para modificar los punteros
    Octree temp_tree = *host_tree;

    // 3. Reservar memoria para cada array en GPU
    float *d_center_x, *d_center_y, *d_center_z, *d_half_size, *d_mass;
    double *d_com_x, *d_com_y, *d_com_z;
    unsigned int (*d_children)[8];
    long *d_star_index;

    size_t size = host_tree->size;
    
    // Verificar que size sea válido
    if (size == 0) {
        printf("Error: tamaño del árbol es 0\n");
        cudaFree(gpu_tree_struct);
        *d_tree = NULL;
        return;
    }

    // Reservar memoria para todos los arrays con verificación de errores
    if ((err = cudaMalloc(&d_center_x, size * sizeof(float))) != cudaSuccess ||
        (err = cudaMalloc(&d_center_y, size * sizeof(float))) != cudaSuccess ||
        (err = cudaMalloc(&d_center_z, size * sizeof(float))) != cudaSuccess ||
        (err = cudaMalloc(&d_half_size, size * sizeof(float))) != cudaSuccess ||
        (err = cudaMalloc(&d_mass, size * sizeof(float))) != cudaSuccess ||
        (err = cudaMalloc(&d_com_x, size * sizeof(double))) != cudaSuccess ||
        (err = cudaMalloc(&d_com_y, size * sizeof(double))) != cudaSuccess ||
        (err = cudaMalloc(&d_com_z, size * sizeof(double))) != cudaSuccess ||
        (err = cudaMalloc(&d_children, size * sizeof(unsigned int[8]))) != cudaSuccess ||
        (err = cudaMalloc(&d_star_index, size * sizeof(long))) != cudaSuccess) {
        
        printf("Error en cudaMalloc para arrays: %s\n", cudaGetErrorString(err));
        
        // Limpiar memoria ya reservada
        cudaFree(d_center_x); cudaFree(d_center_y); cudaFree(d_center_z);
        cudaFree(d_half_size); cudaFree(d_mass); cudaFree(d_com_x);
        cudaFree(d_com_y); cudaFree(d_com_z); cudaFree(d_children);
        cudaFree(d_star_index); cudaFree(gpu_tree_struct);
        *d_tree = NULL;
        return;
    }

    // 4. Copiar datos de los arrays
    cudaMemcpyAsync(d_center_x, host_tree->center_x, size * sizeof(float),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_center_y, host_tree->center_y, size * sizeof(float),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_center_z, host_tree->center_z, size * sizeof(float),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_half_size, host_tree->half_size, size * sizeof(float),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_mass, host_tree->mass, size * sizeof(float),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_com_x, host_tree->com_x, size * sizeof(double),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_com_y, host_tree->com_y, size * sizeof(double),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_com_z, host_tree->com_z, size * sizeof(double),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_children, host_tree->children, size * sizeof(unsigned int[8]),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_star_index, host_tree->star_index, size * sizeof(long),
                    cudaMemcpyHostToDevice, stream);

    // 5. Actualizar los punteros en la estructura temporal
    temp_tree.center_x = d_center_x;
    temp_tree.center_y = d_center_y;
    temp_tree.center_z = d_center_z;
    temp_tree.half_size = d_half_size;
    temp_tree.mass = d_mass;
    temp_tree.com_x = d_com_x;
    temp_tree.com_y = d_com_y;
    temp_tree.com_z = d_com_z;
    temp_tree.children = d_children;
    temp_tree.star_index = d_star_index;

    // 6. Copiar la estructura modificada a GPU
    cudaMemcpyAsync(gpu_tree_struct, &temp_tree, sizeof(Octree),
                    cudaMemcpyHostToDevice, stream);

    // 7. Devolver el puntero a la estructura en GPU
    *d_tree = gpu_tree_struct;
}

void free_tree_gpu(Octree *d_tree) {
    if (d_tree == NULL) return;

    // Primero obtener la estructura desde GPU para liberar los arrays
    Octree temp;
    cudaMemcpy(&temp, d_tree, sizeof(Octree), cudaMemcpyDeviceToHost);

    // Liberar todos los arrays
    cudaFree(temp.center_x);
    cudaFree(temp.center_y);
    cudaFree(temp.center_z);
    cudaFree(temp.half_size);
    cudaFree(temp.mass);
    cudaFree(temp.com_x);
    cudaFree(temp.com_y);
    cudaFree(temp.com_z);
    cudaFree(temp.children);
    cudaFree(temp.star_index);

    // Liberar la estructura principal
    cudaFree(d_tree);
}

__host__ int compute_acceleration_multi_gpu(unsigned int N, int iterations, double *ax, double *ay, double *az,
                                             const unsigned int *offsets, int device_count, Octree **trees,
                                             const Star *estrellas, cudaStream_t *streams) {

    for (int i = 0; i < iterations; i++) {
        #pragma omp parallel for num_threads(device_count)
        for (int dev = 0; dev < device_count; dev++) {
            int octant = i * device_count + dev;
            if (octant >= 8) continue;
            
            cudaError_t err = cudaSetDevice(dev);
            if (err != cudaSuccess) {
                printf("Error setting device %d: %s\n", dev, cudaGetErrorString(err));
                continue;
            }

            long start = offsets[octant];
            long end = (octant == 7) ? N : offsets[octant + 1];
            long count = end - start;

            if (count <= 0) {
                printf("Saltando octante %d: count=%ld\n", octant, count);
                continue;
            }

            // Verificar que el árbol existe
            Octree *tree = trees[octant];
            if (tree == NULL) {
                printf("Error: árbol nulo para octante %d\n", octant);
                continue;
            }
            // Copiar árbol a GPU
            Octree *d_tree;
            copy_tree_to_gpu(&d_tree, tree, streams[dev]);
            if (d_tree == NULL) {
                printf("Error copiando árbol a GPU\n");
                continue;
            }

            // Reservar memoria para coordenadas y aceleraciones
            double *d_x, *d_y, *d_z, *d_ax, *d_ay, *d_az;
            
            if ((err = cudaMalloc(&d_x, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_y, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_z, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_ax, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_ay, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_az, count * sizeof(double))) != cudaSuccess) {
                
                printf("Error reservando memoria coordenadas: %s\n", cudaGetErrorString(err));
                free_tree_gpu(d_tree);
                continue;
            }
            printf("Memoria reservada en GPU it: %d dev: %d count: %ld\n", i, dev, count);
            fflush(stdout);

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

            // CRÍTICO: Sincronizar el stream antes de lanzar el kernel
            cudaStreamSynchronize(streams[dev]);

            printf("Kernel iniciado it:%d dev:%d\n", i, dev); fflush(stdout);

            // Lanzar kernel
            unsigned int grid_size = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
            compute_forces_kernel<<<grid_size, BLOCK_SIZE, 0, streams[dev]>>>(
                d_tree, count, d_x, d_y, d_z, d_ax, d_ay, d_az, 0.2);

            // Verificar errores del kernel
            err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("Error en kernel it:%d dev:%d: %s\n", i, dev,cudaGetErrorString(err));
                fflush(stdout);
                cudaFree(d_x); cudaFree(d_y); cudaFree(d_z);
                cudaFree(d_ax); cudaFree(d_ay); cudaFree(d_az);
                free_tree_gpu(d_tree);
                continue;
            }

            // Sincronizar antes de copiar resultados
            cudaStreamSynchronize(streams[dev]);
            printf("Kernel terminado it:%d dev:%d\n", i, dev); fflush(stdout);

            // Copiar resultados
            cudaMemcpyAsync(&ax[start], d_ax, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);
            cudaMemcpyAsync(&ay[start], d_ay, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);
            cudaMemcpyAsync(&az[start], d_az, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);

            cudaStreamSynchronize(streams[dev]);

            // Liberar memoria
            cudaFree(d_x); cudaFree(d_y); cudaFree(d_z);
            cudaFree(d_ax); cudaFree(d_ay); cudaFree(d_az);
            free_tree_gpu(d_tree);
        }
        
        // Sincronizar todos los dispositivos
        for (int dev = 0; dev < device_count; dev++) {
            cudaSetDevice(dev);
            cudaStreamSynchronize(streams[dev]);
        }
    }
    return 0;
}

// Función principal de simulacion en gpus
extern "C" void simulate_multi_gpu_unified(Star *estrellas,const int steps, const long N, const char *outputfile) {
    struct timeval start, end;
    float cx,cy,cz;
    float hs,min_node_size;
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

        compute_root_bounds(estrellas,&cx,&cy,&cz,&hs,&min_node_size,MIN_SUBDIVISIONS);

        unsigned int offsets[8];
        reorder_stars(estrellas, cx,cy,cz, offsets);
        // Construir árbol
        Octree **octrees = build_tree_gpu(estrellas,cx,cy,cz,hs,min_node_size,offsets);
        printf("Iniciando fase 1\n"); fflush(stdout);
        if (compute_acceleration_multi_gpu(N, iterations, ax, ay, az, offsets, device_count, octrees, estrellas, streams)!=0) {
            printf("Error en fase 1\n");
            return;
        }
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
        for (int i=0;i<8;i++) {
            free_tree(octrees[i]);
        }
        compute_root_bounds(estrellas,&cx,&cy,&cz,&hs,&min_node_size,MIN_SUBDIVISIONS);
        reorder_stars(estrellas, cx,cy,cz, offsets);
        octrees = build_tree_gpu(estrellas,cx,cy,cz,hs,min_node_size,offsets);
        printf("Iniciando fase 2\n"); fflush(stdout);
        if (compute_acceleration_multi_gpu(N, iterations, ax, ay, az, offsets, device_count, octrees, estrellas, streams)!=0) {
            printf("Error en fase 2\n");
            return;
        }
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
    write_chunks(estrellas,"results_cuda",outputfile,25,0);
}