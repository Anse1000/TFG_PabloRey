#include "cuda_functions.h"
#include <cuda_runtime.h>

#define BLOCK_SIZE 256


__global__ void cuda_kernel(const Star *estrellas, double *ax, double *ay, double *az, int N) {
    __shared__ Star sh_stars[BLOCK_SIZE];

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double ax_i = 0.0, ay_i = 0.0, az_i = 0.0;
    double Ci_x = estrellas[i].C[0];
    double Ci_y = estrellas[i].C[1];
    double Ci_z = estrellas[i].C[2];

    for (int tile = 0; tile < (N + BLOCK_SIZE - 1) / BLOCK_SIZE; tile++) {
        int j = tile * BLOCK_SIZE + threadIdx.x;
        if (j < N) {
            sh_stars[threadIdx.x] = estrellas[j];
        }
        __syncthreads();

        #pragma unroll
        for (int k = 0; k < BLOCK_SIZE; k++) {
            int j_global = tile * BLOCK_SIZE + k;
            if (j_global >= N || j_global == i) continue;

            double dx = sh_stars[k].C[0] - Ci_x;
            double dy = sh_stars[k].C[1] - Ci_y;
            double dz = sh_stars[k].C[2] - Ci_z;
            double dist_sq = dx * dx + dy * dy + dz * dz;

            double inv_dist = rsqrt(dist_sq);
            double force = -G * sh_stars[k].mass * inv_dist * inv_dist * inv_dist;

            ax_i += force * dx;
            ay_i += force * dy;
            az_i += force * dz;
        }
    }
    ax[i] = ax_i;
    ay[i] = ay_i;
    az[i] = az_i;
}

extern "C" void compute_aceleration_CUDA(Star *estrellas, double *ax, double *ay, double *az, int N) {
    Star *d_estrellas;
    double *d_ax, *d_ay, *d_az;

    cudaMalloc(&d_estrellas, N * sizeof(Star));
    cudaMalloc(&d_ax, N * sizeof(double));
    cudaMalloc(&d_ay, N * sizeof(double));
    cudaMalloc(&d_az, N * sizeof(double));

    cudaMemcpy(d_estrellas, estrellas, N * sizeof(Star), cudaMemcpyHostToDevice);

    int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    cuda_kernel<<<numBlocks, BLOCK_SIZE>>>(d_estrellas, d_ax, d_ay, d_az, N);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
        return;
    }

    cudaMemcpy(ax, d_ax, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ay, d_ay, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(az, d_az, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_estrellas);
    cudaFree(d_ax);
    cudaFree(d_ay);
    cudaFree(d_az);
}


