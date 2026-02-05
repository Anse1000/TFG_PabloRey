#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

#include "../types.h"

#ifdef __cplusplus
extern "C" {
#endif
    void simulate_multi_gpu_unified(Star *estrellas,int steps,long N, const char *outputfile, float DT);
#ifdef DEBUG_BUILD
    void mem_test_gpu(Star *estrellas);
    void reversibility_test(Star *estrellas);
    void benchmark_tree_construction_GPU(Star *estrellas, int iteraciones);
#endif
#ifdef __cplusplus
}
#endif

#endif //CUDA_FUNCTIONS_H
