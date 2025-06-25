#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H
#include "types.h"

typedef struct {
    double center_x, center_y, center_z;
    double half_size;
    float  mass;
    double com_x, com_y, com_z;
    long   star_index;
} NodoResumen;

#ifdef __cplusplus
extern "C" {
#endif
    void distribute_root(const Octree *tree, NodoResumen **hermanos_gpus);
    void distribute_tree_gpu(const Octree *tree_cpu, int gpu, Octree **tree_gpu_out);
#ifdef __cplusplus
}
#endif
#endif // CUDA_FUNCTIONS_H
