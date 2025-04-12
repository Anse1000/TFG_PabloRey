#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H
#include "types.h"
#ifdef __cplusplus
extern "C" {
#endif
    void compute_aceleration_CUDA(Star *estrellas, double *ax, double *ay, double *az, int N);
#ifdef __cplusplus
}
#endif


#endif // CUDA_FUNCTIONS_H
