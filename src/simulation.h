#ifndef SIMULATION_H
#define SIMULATION_H
#include "types.h"
#ifdef AVX_512
#include <immintrin.h>
#endif
#ifdef CUDA
#include "cuda_functions.h"
#endif

typedef struct OctreeNode {
    double cx, cy, cz;      // centro del cubo
    float half_size;       // mitad del tamaño del cubo

    float mass;            // masa total dentro del nodo
    double com_x, com_y, com_z; // centro de masa

    short is_leaf;
    long star_index;

    struct OctreeNode *children[8];  // hijos (pueden ser NULL)
} OctreeNode;


void simulate(Star *estrellas, int N);
void test_simulation(Star *estrellas);
#endif //SIMULATION_H
