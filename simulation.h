#ifndef SIMULATION_H
#define SIMULATION_H
#include "types.h"
#ifdef AVX_512
#include <immintrin.h>
#endif
void simulate(Star *estrellas, const int N);
#endif //SIMULATION_H
