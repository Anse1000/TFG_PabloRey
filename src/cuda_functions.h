#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H
#include "types.h"
#ifdef __cplusplus
extern "C" {
#endif
    void compute_aceleration_CUDA(Star *estrellas, float *ax, float *ay, float *az, int N);
#ifdef __cplusplus
}
#endif


#endif // CUDA_FUNCTIONS_H
