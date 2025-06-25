#include <cuda_runtime.h>
#include "cuda_functions.h"

void print_gpu_memory_info(int device) {
    cudaSetDevice(device);

    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);

    double free_MB = free_mem / (1024.0 * 1024.0);
    double total_MB = total_mem / (1024.0 * 1024.0);
    double used_MB = total_MB - free_MB;

    printf("GPU %d: Used %.2f MB / %.2f MB (%.2f%% used)\n",
           device, used_MB, total_MB, 100.0 * used_MB / total_MB);
}

extern "C" void distribute_root(const Octree *tree, NodoResumen **hermanos_gpus) {
    long *hijos = tree->children[0];  // Hijos del nodo raíz
    if (!hijos) {
        fprintf(stderr, "ERROR: tree->children[0] es NULL\n");
        exit(1);
    }
    NodoResumen resumenes[8];

    // Construimos los 8 resúmenes en CPU
    for (int i = 0; i < 8; ++i) {
        long idx = hijos[i];
        if (idx == -1) {
            fprintf(stderr, "Error: hijo raíz vacío en GPU %d\n", i);
            exit(1);
        }
        resumenes[i].center_x   = tree->center_x[idx];
        resumenes[i].center_y   = tree->center_y[idx];
        resumenes[i].center_z   = tree->center_z[idx];
        resumenes[i].half_size  = tree->half_size[idx];
        resumenes[i].mass       = tree->mass[idx];
        resumenes[i].com_x      = tree->com_x[idx];
        resumenes[i].com_y      = tree->com_y[idx];
        resumenes[i].com_z      = tree->com_z[idx];
        resumenes[i].star_index = tree->star_index[idx];
    }

    // Distribuir a cada GPU los 7 hermanos
    for (int i = 0; i < 8; ++i) {
        cudaSetDevice(i);

        NodoResumen hermanos[7];
        int count = 0;
        for (int j = 0; j < 8; ++j) {
            if (j == i) continue;
            hermanos[count++] = resumenes[j];
        }

        NodoResumen *d_hermanos;
        cudaMalloc(&d_hermanos, 7 * sizeof(NodoResumen));
        cudaMemcpy(d_hermanos, hermanos, 7 * sizeof(NodoResumen), cudaMemcpyHostToDevice);

        hermanos_gpus[i] = d_hermanos;
    }
    printf("Hermanos distribuidos a GPU\n"); fflush(stdout);
}

extern "C" void distribute_tree_gpu(const Octree *tree_cpu, int gpu, Octree **tree_gpu_out) {
    long root = tree_cpu->children[0][gpu];
    if (root == -1) {
        fprintf(stderr, "No root node for GPU %d\n", gpu);
        exit(1);
    }

    cudaSetDevice(gpu);
    size_t size = tree_cpu->size/4;

    // Paso 1: Recorrer subárbol (desde CPU)
    long *queue = (long*)malloc(size * sizeof(long));  // Queue: global_idx (nodo
    long *global_to_local= (long*)malloc(size * sizeof(long));  // Map: global_idx → local_idx
    memset(global_to_local, -1, size * sizeof(long));

    size_t count = 0;
    queue[count] = root;
    global_to_local[root] = 0;
    count++;

    for (size_t i = 0; i < count; ++i) {
        long curr = queue[i];
        for (int j = 0; j < 8; ++j) {
            long child = tree_cpu->children[curr][j];
            if (child != -1 && global_to_local[child] == -1) {
                global_to_local[child] = count;
                queue[count++] = child;
            }
        }
    }

    size_t total = count;

    // Paso 2: Reservar memoria en GPU directamente
    Octree *tree_gpu;
    cudaMalloc(&tree_gpu, sizeof(Octree));

    #define GPU_ALLOC(field, type) \
        type *d_##field; \
        cudaMalloc(&d_##field, sizeof(type) * total); \
        cudaMemcpy(&(tree_gpu->field), &d_##field, sizeof(type*), cudaMemcpyHostToDevice);

    GPU_ALLOC(center_x, double);
    GPU_ALLOC(center_y, double);
    GPU_ALLOC(center_z, double);
    GPU_ALLOC(half_size, double);
    GPU_ALLOC(mass, float);
    GPU_ALLOC(com_x, double);
    GPU_ALLOC(com_y, double);
    GPU_ALLOC(com_z, double);
    GPU_ALLOC(star_index, long);

    long (*d_children)[8];
    cudaMalloc(&d_children, sizeof(long[8]) * total);
    cudaMemcpy(&tree_gpu->children, &d_children, sizeof(long(*)[8]), cudaMemcpyHostToDevice);

    cudaMemcpy(&tree_gpu->size, &total, sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(&tree_gpu->capacity, &total, sizeof(size_t), cudaMemcpyHostToDevice);

    // Paso 3: Copiar datos nodo a nodo desde CPU → GPU directo
    for (size_t i = 0; i < total; ++i) {
        long g = queue[i];

        cudaMemcpy(&tree_gpu->center_x[i], &tree_cpu->center_x[g], sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->center_y[i], &tree_cpu->center_y[g], sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->center_z[i], &tree_cpu->center_z[g], sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->half_size[i], &tree_cpu->half_size[g], sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->mass[i],       &tree_cpu->mass[g],      sizeof(float),  cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->com_x[i],      &tree_cpu->com_x[g],     sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->com_y[i],      &tree_cpu->com_y[g],     sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->com_z[i],      &tree_cpu->com_z[g],     sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&tree_gpu->star_index[i], &tree_cpu->star_index[g],sizeof(long),   cudaMemcpyHostToDevice);

        // Procesar hijos
        long local_children[8];
        for (int j = 0; j < 8; ++j) {
            long c = tree_cpu->children[g][j];
            local_children[j] = c == -1 ? -1 : global_to_local[c];
        }
        cudaMemcpy(&tree_gpu->children[i], local_children, sizeof(long[8]), cudaMemcpyHostToDevice);
    }

    *tree_gpu_out = tree_gpu;
    print_gpu_memory_info(gpu);

}






