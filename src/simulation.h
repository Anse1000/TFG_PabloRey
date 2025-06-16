#ifndef SIMULATION_H
#define SIMULATION_H
#include "types.h"
#ifdef AVX_512
#include <immintrin.h>
#endif
#ifdef CUDA
#include "cuda_functions.h"
#endif

void simulate(Star *estrellas, int N);
void test_simulation(Star *estrellas);
void predict_tree(Star *estrellas);
#endif //SIMULATION_H
