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
#endif //SIMULATION_H
