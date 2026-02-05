#ifndef TFG_PRA_OCTREE_GPU_H
#define TFG_PRA_OCTREE_GPU_H

#include <stdio.h>
#include <stdint.h>
#include "../types.h"

#define INVALID_INDEX UINT32_MAX

typedef struct {
    float center_x, center_y, center_z;
    float half_size;
    float mass;
    double com_x, com_y, com_z;
    unsigned int children[8];
    long star_index;
    unsigned int next;
} OctreeNode;

typedef struct {
    OctreeNode *nodes;
    size_t size;
    size_t capacity;
} OctreeOctant;

typedef struct {
    double com_x, com_y, com_z;
    float mass;
} FrontierNode;

typedef struct {
    FrontierNode *nodes;
    size_t size;
    size_t capacity;
} FrontierGPU;

typedef struct {
    OctreeOctant *octants[8];
    FrontierGPU *frontier[8];
} OctreeGPU;

#ifdef __cplusplus
extern "C" {
#endif
    OctreeGPU *build_tree_GPU(Star *stars, float cx, float cy, float cz, float hs, float min_node_size,const unsigned int *offsets);
    void free_tree_gpu(OctreeGPU *tree);
#ifdef __cplusplus
}
#endif
#endif //TFG_PRA_OCTREE_GPU_H