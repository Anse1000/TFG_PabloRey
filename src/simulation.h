#ifndef SIMULATION_H
#define SIMULATION_H
#include "types.h"
#ifdef CUDA
#include "cuda_functions.h"
#endif

#define INVALID_INDEX UINT32_MAX

void simulate(Star *estrellas, long N, const char* outputfile);
void test_simulation(Star *estrellas);

#endif //SIMULATION_H
