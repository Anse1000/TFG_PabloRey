#ifndef OCTREE_H
#define OCTREE_H

#include <stdint.h>
#include <stdio.h>
#include "types.h"

#define INVALID_INDEX UINT32_MAX

typedef struct {
    // Bounding box (centro y tamaño)
    float *center_x, *center_y, *center_z;
    float *half_size;

    // Agregado de masa
    float *mass;
    double *com_x, *com_y, *com_z; // centro de masa

    // Hijos (índices, -1 si no existe). 8 hijos por nodo.
    unsigned int (*children)[8]; // tamaño = capacity

    // Índice de estrella si hoja, -1 si nodo interno
    long *star_index;

    size_t size; // nodos usados
    size_t capacity; // capacidad total
} Octree;
#ifdef __cplusplus
extern "C" {
#endif
Octree *build_tree(Star *stars);
#ifdef CUDA
Octree **build_tree_gpu(Star *stars, float cx, float cy, float cz, float hs, float min_node_size,
                        const unsigned int *offsets);
#endif
#ifdef DEBUG_BUILD
void test_tree(Star *stars);
#endif
#ifdef __cplusplus
}
#endif
#endif //OCTREE_H
