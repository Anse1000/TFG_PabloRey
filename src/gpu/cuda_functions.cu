#include "cuda_functions.cuh"
#include "../aux_fun.h"
#include <sys/time.h>
#include "octree_gpu.h"
#include "../aux_fun.h"
#include "../file_handler.h"

#define BLOCK_SIZE 256

__constant__ double c_dt;
__constant__ double c_dt2;
__constant__ float c_theta2;
__constant__ size_t c_star_count;
__constant__ long c_base_index;

// --- Función para aceleración del halo NFW
__device__ __forceinline__ double halo_accel_gpu(double r, double *ax, double *ay, double *az,
                                                 double dx, double dy, double dz) {
    double x = r / rs;
    double f = log1p(x) - x / (1.0 + x);
    double Menc = M200 * f / 1.488804364; //(log(1.0 + 10.0) - 10.0 / 11.0)
    double inv_r = 1.0 / r;
    double inv_r3 = inv_r * inv_r * inv_r;
    double acc = -G * Menc * inv_r3;

    // FMA para precisión y rendimiento
    *ax = fma(acc, dx, *ax);
    *ay = fma(acc, dy, *ay);
    *az = fma(acc, dz, *az);

    return acc;
}

// --- Función para aceleración del bulbo Hernquist
__device__ __forceinline__ double bulge_accel_gpu(double r, double *ax, double *ay, double *az,
                                                  double dx, double dy, double dz) {
    // -G M / (r (r + A)^2)
    double rpA = r + A;
    double inv_r = 1.0 / r;
    double inv_rpA2 = 1.0 / (rpA * rpA);
    double acc = -G * MBULGE * inv_r * inv_rpA2;


    *ax = fma(acc, dx, *ax);
    *ay = fma(acc, dy, *ay);
    *az = fma(acc, dz, *az);

    return acc;
}

// --- Función wrapper para aceleración analítica total
__device__ __forceinline__ void analytic_accel_gpu(double x, double y, double z,
                                                   double *ax, double *ay, double *az) {
    double dx = x;
    double dy = y;
    double dz = z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    if (r > 0) {
        halo_accel_gpu(r, ax, ay, az, dx, dy, dz);
        bulge_accel_gpu(r, ax, ay, az, dx, dy, dz);
    }
}

__global__ void compute_forces_kernel(const OctreeOctant *__restrict__ tree,
                                      const FrontierGPU *__restrict__ frontier,
                                      double *__restrict__ Cx, double *__restrict__ Cy, double *__restrict__ Cz,
                                      double *__restrict__ Vx, double *__restrict__ Vy, double *__restrict__ Vz,
                                      const bool do_drift) {
    unsigned int star_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (star_idx >= c_star_count) return;
    long global_idx = c_base_index + star_idx;

    // Inicializar aceleraciones
    double acc_x = 0.0, acc_y = 0.0, acc_z = 0.0;

    // --- Aceleración analítica
    analytic_accel_gpu(Cx[star_idx], Cy[star_idx], Cz[star_idx],
                       &acc_x, &acc_y, &acc_z);

    double dx, dy, dz, dist_sq, inv_dist, inv_dist3, force;

    // --- Traversal stackless del árbol local usando ropes
    unsigned int node_idx = 0; // raíz del octante
    while (node_idx != INVALID_INDEX) {
        OctreeNode node = tree->nodes[node_idx];
        // Saltar si es la misma estrella
        if (node.star_index != global_idx) {
            dx = node.com_x - Cx[star_idx];
            dy = node.com_y - Cy[star_idx];
            dz = node.com_z - Cz[star_idx];
            dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;

            float s = 2.0F * node.half_size;

            if (s * s < c_theta2 * dist_sq || node.star_index >= 0) {
                inv_dist = rsqrt(dist_sq);
                inv_dist3 = inv_dist * inv_dist * inv_dist;
                force = -G * node.mass * inv_dist3;

                acc_x = fma(force, dx, acc_x);
                acc_y = fma(force, dy, acc_y);
                acc_z = fma(force, dz, acc_z);

                // Saltar usando rope
                node_idx = node.next;
                continue;
            }
        }

        // Bajar al primer hijo si existe
        unsigned int child = INVALID_INDEX;
        for (int i = 0; i < 8; i++) {
            if (tree->nodes[node_idx].children[i] != INVALID_INDEX) {
                child = tree->nodes[node_idx].children[i];
                break;
            }
        }

        if (child != INVALID_INDEX) {
            node_idx = child;
        } else {
            // Saltar usando rope
            node_idx = tree->nodes[node_idx].next;
        }
    }

    // --- Recorrer frontera lineal (ya filtrada en CPU)
    for (size_t i = 0; i < frontier->size; i++) {
        dx = frontier->nodes[i].com_x - Cx[star_idx];
        dy = frontier->nodes[i].com_y - Cy[star_idx];
        dz = frontier->nodes[i].com_z - Cz[star_idx];
        dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;

        inv_dist = rsqrt(dist_sq);
        inv_dist3 = inv_dist * inv_dist * inv_dist;
        force = -G * frontier->nodes[i].mass * inv_dist3;

        acc_x = fma(force, dx, acc_x);
        acc_y = fma(force, dy, acc_y);
        acc_z = fma(force, dz, acc_z);
    }

    // --- Actualizar velocidades
    Vx[star_idx] = fma(c_dt2, acc_x, Vx[star_idx]);
    Vy[star_idx] = fma(c_dt2, acc_y, Vy[star_idx]);
    Vz[star_idx] = fma(c_dt2, acc_z, Vz[star_idx]);

    // --- Drift opcional
    if (do_drift) {
        Cx[star_idx] = fma(c_dt, Vx[star_idx], Cx[star_idx]);
        Cy[star_idx] = fma(c_dt, Vy[star_idx], Cy[star_idx]);
        Cz[star_idx] = fma(c_dt, Vz[star_idx], Cz[star_idx]);
    }
}


static cudaError_t copy_tree_to_gpu(OctreeOctant **d_oct_out, const OctreeOctant *h_oct, FrontierGPU **d_front_out,
                                    FrontierGPU *h_front, cudaStream_t stream) {
    if (!h_oct || !d_oct_out || !h_front || !d_front_out) return cudaErrorInvalidValue;

    cudaError_t err;

    // Copiar nodos del octante
    OctreeNode *d_nodes = nullptr;
    if ((err = cudaMalloc(&d_nodes, h_oct->size * sizeof(OctreeNode))) != cudaSuccess) return err;
    if ((err = cudaMemcpyAsync(d_nodes, h_oct->nodes, h_oct->size * sizeof(OctreeNode),
                               cudaMemcpyHostToDevice, stream)) != cudaSuccess)
        return err;

    // Construir estructura OctreeOctant en host con puntero de dispositivo
    OctreeOctant h_oct_dev = {};
    h_oct_dev.nodes = d_nodes;
    h_oct_dev.size = h_oct->size;
    h_oct_dev.capacity = h_oct->size;

    // Reservar y copiar estructura OctreeOctant al dispositivo
    OctreeOctant *d_oct = nullptr;
    if ((err = cudaMalloc(&d_oct, sizeof(OctreeOctant))) != cudaSuccess) return err;
    if ((err = cudaMemcpyAsync(d_oct, &h_oct_dev, sizeof(OctreeOctant),
                               cudaMemcpyHostToDevice, stream)) != cudaSuccess)
        return err;

    // Copiar frontera
    FrontierNode *d_fnodes = nullptr;
    if ((err = cudaMalloc(&d_fnodes, h_front->size * sizeof(FrontierNode))) != cudaSuccess) return err;
    if ((err = cudaMemcpyAsync(d_fnodes, h_front->nodes, h_front->size * sizeof(FrontierNode),
                               cudaMemcpyHostToDevice, stream)) != cudaSuccess)
        return err;

    FrontierGPU h_front_dev = {};
    h_front_dev.nodes = d_fnodes;
    h_front_dev.size = h_front->size;
    h_front_dev.capacity = h_front->size;

    FrontierGPU *d_front = nullptr;
    if ((err = cudaMalloc(&d_front, sizeof(FrontierGPU))) != cudaSuccess) return err;
    if ((err = cudaMemcpyAsync(d_front, &h_front_dev, sizeof(FrontierGPU),
                               cudaMemcpyHostToDevice, stream)) != cudaSuccess)
        return err;

    *d_oct_out = d_oct;
    *d_front_out = d_front;
    return cudaSuccess;
}

// Liberación de memoria de dispositivo para octante y frontera
static void free_octant_and_frontier_gpu(OctreeOctant *d_oct, FrontierGPU *d_front) {
    if (d_oct) {
        // Recuperar la estructura a host para liberar los buffers internos
        OctreeOctant h_oct_dev;
        if (cudaMemcpy(&h_oct_dev, d_oct, sizeof(OctreeOctant), cudaMemcpyDeviceToHost) == cudaSuccess) {
            if (h_oct_dev.nodes) cudaFree(h_oct_dev.nodes);
        }
        cudaFree(d_oct);
    }
    if (d_front) {
        FrontierGPU h_front_dev;
        if (cudaMemcpy(&h_front_dev, d_front, sizeof(FrontierGPU), cudaMemcpyDeviceToHost) == cudaSuccess) {
            if (h_front_dev.nodes) cudaFree(h_front_dev.nodes);
        }
        cudaFree(d_front);
    }
}

__host__ int compute_halfstep(unsigned int N, int iterations, const unsigned int *offsets, int device_count,
                              OctreeGPU *tree, const Star *estrellas, cudaStream_t *streams, const double DT,
                              int drift) {
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
            OctreeOctant *octtree = tree->octants[octant];
            FrontierGPU *frontier = tree->frontier[octant];
            if (octtree == NULL || frontier == NULL) {
                printf("Error: árbol nulo para octante %d\n", octant);
                continue;
            }
            // Copiar arbol a GPU
            OctreeOctant *d_tree = NULL;
            FrontierGPU *d_front = NULL;
            if ((err = copy_tree_to_gpu(&d_tree, octtree, &d_front, frontier, streams[dev]))) {
                printf("Error copiando árbol/frontera a GPU: %s\n", cudaGetErrorString(err));
                continue;
            }

            // Reservar memoria para coordenadas, velocidades y constantes
            double *d_cx, *d_cy, *d_cz;
            double *d_vx, *d_vy, *d_vz;

            double DT2 = DT * 0.5;
            float theta = THETA;
            float theta2 = theta * theta;

            cudaMemcpyToSymbol(c_dt, &DT, sizeof(double));
            cudaMemcpyToSymbol(c_dt2, &DT2, sizeof(double));
            cudaMemcpyToSymbol(c_theta2, &theta2, sizeof(float));
            cudaMemcpyToSymbol(c_star_count, &count, sizeof(size_t));
            cudaMemcpyToSymbol(c_base_index, &start, sizeof(long));


            if ((err = cudaMalloc(&d_cx, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_cy, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_cz, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_vx, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_vy, count * sizeof(double))) != cudaSuccess ||
                (err = cudaMalloc(&d_vz, count * sizeof(double))) != cudaSuccess) {
                printf("Error reservando memoria en gpu: %s\n", cudaGetErrorString(err));
                free_octant_and_frontier_gpu(d_tree, d_front);
                continue;
            }
#ifdef DEBUG_BUILD
            printf("Memoria reservada en GPU it: %d dev: %d count: %ld\n", i, dev, count);
            fflush(stdout);
#endif

            // Copiar datos a GPU
            cudaMemcpyAsync(d_cx, &estrellas->Cx[start], count * sizeof(double),
                            cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_cy, &estrellas->Cy[start], count * sizeof(double),
                            cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_cz, &estrellas->Cz[start], count * sizeof(double),
                            cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_vx, &estrellas->Vx[start], count * sizeof(double),
                            cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_vy, &estrellas->Vy[start], count * sizeof(double),
                            cudaMemcpyHostToDevice, streams[dev]);
            cudaMemcpyAsync(d_vz, &estrellas->Vz[start], count * sizeof(double),
                            cudaMemcpyHostToDevice, streams[dev]);
#ifdef DEBUG_BUILD
            printf("Kernel iniciado it:%d dev:%d\n", i, dev);
            fflush(stdout);
#endif
            // Lanzar kernels
            unsigned int grid_size = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;
            compute_forces_kernel<<<grid_size, BLOCK_SIZE, 0, streams[dev]>>>(
                d_tree, d_front, d_cx, d_cy, d_cz, d_vx, d_vy, d_vz, drift);
            // Verificar errores del kernel
            err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("Error en kernel it:%d dev:%d: %s\n", i, dev, cudaGetErrorString(err));
                fflush(stdout);
                cudaFree(d_cx);
                cudaFree(d_cy);
                cudaFree(d_cz);
                cudaFree(d_vx);
                cudaFree(d_vy);
                cudaFree(d_vz);
                free_octant_and_frontier_gpu(d_tree, d_front);
                continue;
            }

            cudaStreamSynchronize(streams[dev]);
#ifdef DEBUG_BUILD
            printf("Kernel terminado it:%d dev:%d\n", i, dev);
            fflush(stdout);
#endif
            // Copiar datos a CPU
            if (drift) {
                cudaMemcpyAsync(&estrellas->Cx[start], d_cx, count * sizeof(double),
                                cudaMemcpyDeviceToHost, streams[dev]);
                cudaMemcpyAsync(&estrellas->Cy[start], d_cy, count * sizeof(double),
                                cudaMemcpyDeviceToHost, streams[dev]);
                cudaMemcpyAsync(&estrellas->Cz[start], d_cz, count * sizeof(double),
                                cudaMemcpyDeviceToHost, streams[dev]);
            }
            cudaMemcpyAsync(&estrellas->Vx[start], d_vx, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);
            cudaMemcpyAsync(&estrellas->Vy[start], d_vy, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);
            cudaMemcpyAsync(&estrellas->Vz[start], d_vz, count * sizeof(double),
                            cudaMemcpyDeviceToHost, streams[dev]);

            cudaStreamSynchronize(streams[dev]);
#ifdef DEBUG_BUILD
            printf("Copia de datos terminada it:%d dev:%d\n", i, dev);
            fflush(stdout);
#endif

            // Liberar memoria
            cudaFree(d_cx);
            cudaFree(d_cy);
            cudaFree(d_cz);
            cudaFree(d_vx);
            cudaFree(d_vy);
            cudaFree(d_vz);
            free_octant_and_frontier_gpu(d_tree, d_front);
#ifdef DEBUG_BUILD
            printf("Memoria liberada it:%d dev:%d\n", i, dev);
            fflush(stdout);
#endif
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
extern "C" void simulate_multi_gpu_unified(Star *estrellas, const int steps, const long N, const char *outputfile,
                                           const float DT) {
    struct timeval start, end;
    float cx, cy, cz;
    float hs, min_node_size;
    gettimeofday(&start, NULL);
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        printf("No hay dispositivos disponibles\n");
        exit(1);
    }
    // Crear streams para cada GPU
    auto *streams = static_cast<cudaStream_t *>(malloc(device_count * sizeof(cudaStream_t)));
    for (int i = 0; i < device_count; i++) {
        cudaSetDevice(i);
        cudaStreamCreate(&streams[i]);
    }
    const int iterations = (8 + device_count - 1) / device_count;
    printf("******************************************************\n");
    printf("Iniciando simulacion con %d GPUs\n", device_count);
    printf("Simulando %d pasos de %.0f años (Total: %.0f años)\n", steps, DT * 1000000, DT * steps * 1000000);
    printf("******************************************************\n");
    fflush(stdout);
    compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size,MIN_SUBDIVISIONS);
    unsigned int offsets[8];
    reorder_stars(estrellas, cx, cy, cz, offsets);
    write_results_hdf5(estrellas, outputfile, "cuda_results", -1);
    // Construir árbol
    OctreeGPU *tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
    for (int step = 0; step < steps; step++) {
        struct timeval step_start, step_end;
        printf("****************** Iniciando paso %d ******************\n", step + 1);
        printf("\t  Iniciando fase 1: HalfKick-Drift \n");
        fflush(stdout);
        gettimeofday(&step_start, NULL);

        if (compute_halfstep(N, iterations, offsets, device_count, tree, estrellas,
                             streams, DT, 1) != 0) {
            printf("Error en fase 1\n");
            return;
        }
        free_tree_gpu(tree);
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size,MIN_SUBDIVISIONS);
        reorder_stars(estrellas, cx, cy, cz, offsets);
        tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
        printf("\t  Iniciando fase 2: HalfKick\n");
        fflush(stdout);
        if (compute_halfstep(N, iterations, offsets, device_count, tree, estrellas,
                             streams, DT, 0) != 0) {
            printf("Error en fase 2\n");
            return;
        }
        write_results_hdf5(estrellas, outputfile, "cuda_results", step);
        gettimeofday(&step_end, NULL);
        double step_seconds = get_seconds(step_start, step_end);
        printf("************ Paso %d finalizado en %6.0f segundos ************\n", step + 1, step_seconds);
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
    printf("Simulacion de %ld estrellas completada en %02d:%02d:%05.2f (hh:mm:ss)\n",
           N, hours, minutes, remaining_seconds);
    printf("Resultados guardados en %s\n", outputfile);
    fflush(stdout);
}

#ifdef DEBUG_BUILD
#include <algorithm>

// Estructura auxiliar para no mover todos los datos pesados durante el sort
struct SortPair {
    unsigned long id;
    long original_idx;
};

extern "C" void sort_stars_by_id(Star *estrellas) {
    const long N = (long)estrellas->size;

    // 1. Crear array de pares (ID, índice actual)
    auto *pairs = (SortPair*)malloc(N * sizeof(SortPair));

    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        pairs[i].id = estrellas->id[i];
        pairs[i].original_idx = i;
    }

    // 2. Ordenar los pares por ID
    // Usamos std::sort (que es muy rápido) con una lambda.
    // Nota: Aunque std::sort no es paralelo de serie, es preferible a un quicksort manual.
    // Si tu compilador soporta C++17 y tienes TBB, podrías usar std::sort(std::execution::par, ...)
    std::sort(pairs, pairs + N, [](const SortPair& a, const SortPair& b) {
        return a.id < b.id;
    });

    // 3. Reubicar los datos usando un buffer temporal por array para evitar colisiones
    // Solo necesitamos un buffer temporal por cada tipo de dato (double/unsigned long)
    auto *tmp_double = (double*)malloc(N * sizeof(double));
    auto *tmp_id = (unsigned long*)malloc(N * sizeof(unsigned long));

    auto apply_permutation_double = [&](double* target_array) {
        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            tmp_double[i] = target_array[pairs[i].original_idx];
        }
        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            target_array[i] = tmp_double[i];
        }
    };

    // Aplicar a todos los arrays de posición y velocidad
    apply_permutation_double(estrellas->Cx);
    apply_permutation_double(estrellas->Cy);
    apply_permutation_double(estrellas->Cz);
    apply_permutation_double(estrellas->Vx);
    apply_permutation_double(estrellas->Vy);
    apply_permutation_double(estrellas->Vz);

    // Aplicar al array de IDs
    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        tmp_id[i] = estrellas->id[pairs[i].original_idx];
    }
    #pragma omp parallel for
    for (long i = 0; i < N; i++) {
        estrellas->id[i] = tmp_id[i];
    }
    printf("Estrellas ordenadas por ID\n"); fflush(stdout);
    // Limpiar
    free(pairs);
    free(tmp_double);
    free(tmp_id);
}
extern "C" void mem_test_gpu(Star *estrellas) {
    float cx, cy, cz;
    float hs, min_node_size;
    compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size,MIN_SUBDIVISIONS);
    unsigned int offsets[8];
    reorder_stars(estrellas, cx, cy, cz, offsets);
    // Construir árbol
    OctreeGPU *tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
    for (int i = 0; i < 8; i++) {
        unsigned long start_idx = offsets[i];
        unsigned long end_idx = (i == 7) ? estrellas->size : offsets[i + 1];
        unsigned long count = end_idx - start_idx;
        size_t memory_tree = tree->octants[i]->capacity * (
                                 sizeof(double) * 3 +
                                 sizeof(double) +
                                 sizeof(float) * 4 +
                                 sizeof(unsigned int[8]) +
                                 sizeof(long)
                             ) / 1024 / 1024;
        size_t memory_stars = count * sizeof(double) * 6 / 1024 / 1024;
        size_t memory_frontier = tree->frontier[i]->capacity * sizeof(FrontierNode) / 1024 / 1024;
        printf("Frontier nodes %lu\n", tree->frontier[i]->size);
        printf(
            "Octante %d: %ld estrellas, Memoria estimada: %lu MB subarbol + %lu MB estrellas + %lu MB frontera = %lu MB TOTAL\n",
            i, count, memory_tree, memory_stars, memory_frontier, memory_tree + memory_stars + memory_frontier);
        fflush(stdout);
    }
}

extern "C" void reversibility_test(Star *estrellas) {
    printf("\n========== INICIANDO TEST DE REVERSIBILIDAD TEMPORAL (LEAPFROG KDK) ==========\n");
    const int steps = 10;  // 50 adelante, 50 atrás
    const float DT = 1.0F;
    const long N = (long)estrellas->size;
    int device_count=0;
    cudaGetDeviceCount(&device_count);

    // 1. Backup del estado inicial (t = 0)
    auto *orig_x = (double*)malloc(N * sizeof(double));
    auto *orig_y = (double*)malloc(N * sizeof(double));
    auto *orig_z = (double*)malloc(N * sizeof(double));
    double v_sum = 0;
    sort_stars_by_id(estrellas);
    #pragma omp parallel for reduction(+:v_sum)
    for(long i=0; i<N; i++) {
        orig_x[i] = estrellas->Cx[i];
        orig_y[i] = estrellas->Cy[i];
        orig_z[i] = estrellas->Cz[i];
        v_sum += sqrt(estrellas->Vx[i]*estrellas->Vx[i] + estrellas->Vy[i]*estrellas->Vy[i] + estrellas->Vz[i]*estrellas->Vz[i]);
    }
    double v_avg = v_sum / N;

    if (device_count == 0) {
        printf("No hay dispositivos disponibles\n");
        exit(1);
    }

    auto *streams = (cudaStream_t *)malloc(device_count * sizeof(cudaStream_t));
    for (int i = 0; i < device_count; i++) { cudaSetDevice(i); cudaStreamCreate(&streams[i]); }
    const int iterations = (8 + device_count - 1) / device_count;

    printf("Parámetros: %d pasos | DT = %.2f | Velocidad media: %.4f kpc/Myr\n", steps, DT, v_avg);
    printf("Desplazamiento total esperado por estrella: ~%.4f kpc\n", v_avg * steps * DT);

    // 2. CICLO ADELANTE
    printf("Evolucionando hacia adelante...\n");
    for (int step = 0; step < steps; step++) {
        float cx, cy, cz, hs, min_node_size; unsigned int offsets[8];
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
        reorder_stars(estrellas, cx, cy, cz, offsets);
        OctreeGPU *tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
        compute_halfstep(N, iterations, offsets, device_count, tree, estrellas, streams, DT, 1); // Kick+Drift
        free_tree_gpu(tree);
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
        reorder_stars(estrellas, cx, cy, cz, offsets);
        tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
        compute_halfstep(N, iterations, offsets, device_count, tree, estrellas, streams, DT, 0); // Kick
        free_tree_gpu(tree);
    }

    // 3. CICLO ATRÁS (DT negativo)
    printf("Evolucionando hacia atrás (DT = %.2f)...\n", -DT);
    for (int step = 0; step < steps; step++) {
        float cx, cy, cz, hs, min_node_size; unsigned int offsets[8];
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
        reorder_stars(estrellas, cx, cy, cz, offsets);
        OctreeGPU *tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
        compute_halfstep(N, iterations, offsets, device_count, tree, estrellas, streams, -DT, 1);
        free_tree_gpu(tree);
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
        reorder_stars(estrellas, cx, cy, cz, offsets);
        tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size, offsets);
        compute_halfstep(N, iterations, offsets, device_count, tree, estrellas, streams, -DT, 0);
        free_tree_gpu(tree);
    }
    sort_stars_by_id(estrellas);
    // 4. Análisis de resultados
    double max_err = 0, avg_err = 0;
    #pragma omp parallel for reduction(+:avg_err) reduction(max:max_err)
    for (long i = 0; i < N; i++) {
        double dx = estrellas->Cx[i] - orig_x[i];
        double dy = estrellas->Cy[i] - orig_y[i];
        double dz = estrellas->Cz[i] - orig_z[i];
        double err = sqrt(dx*dx + dy*dy + dz*dz);
        if (err > max_err) max_err = err;
        avg_err += err;
    }
    avg_err /= N;

    printf("\n--- RESULTADOS DEL TEST DE REVERSIBILIDAD ---\n");
    printf("Error máximo de posición:  %.12le kpc\n", max_err);
    printf("Error promedio de posición: %.12le kpc\n", avg_err);
    printf("Deriva por paso (promedio): %.12le kpc/step\n", avg_err / (steps * 2));

    // Comparación con el movimiento real para contexto
    double precision_relativa = avg_err / (v_avg * steps * DT);
    printf("Error relativo al movimiento total: %.4e (Ideal < 1e-10)\n", precision_relativa);
    printf("---------------------------------------------\n\n");

    free(orig_x); free(orig_y); free(orig_z);
    for (int i = 0; i < device_count; i++) { cudaSetDevice(i); cudaStreamDestroy(streams[i]); }
    free(streams);
}
#endif
