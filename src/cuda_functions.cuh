#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif
    void simulate_multi_gpu_unified(Star *estrellas,long N, const char *outputfile);
#ifdef __cplusplus
}
#endif

#endif //CUDA_FUNCTIONS_H
