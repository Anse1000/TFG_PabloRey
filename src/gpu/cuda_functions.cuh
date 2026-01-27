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
#endif
#ifdef __cplusplus
}
#endif

#endif //CUDA_FUNCTIONS_H
