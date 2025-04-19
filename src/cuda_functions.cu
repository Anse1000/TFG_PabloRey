#include "cuda_functions.h"
#include <cuda_runtime.h>

#define BLOCK_SIZE 256


__global__ void cuda_kernel(double *Cx,double *Cy,double *Cz,double *mass, double *ax, double *ay, double *az, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N) {
        ax[i] = 0.0;
        ay[i] = 0.0;
        az[i] = 0.0;

        for (int j = 0; j < N; j++) {
            if (i != j) {
                // Calcular distancia a la estrella
                double dx = Cx[i] - Cx[j];
                double dy = Cy[i] - Cy[j];
                double dz = Cz[i] - Cz[j];
                double dist_sq = dx * dx + dy * dy + dz * dz;
                double inv_dist = rsqrt(dist_sq);

                // Calcular fuerza aplicada a la estrella
                double force = -G * mass[j] * inv_dist * inv_dist * inv_dist;
                ax[i] fma(force, dx, ax[i]);
                ay[i] fma(force, dy, ay[i]);
                az[i] fma(force, dz, az[i]);
            }
        }
    }
}


extern "C" void compute_aceleration_CUDA(Star *stars, double *ax, double *ay, double *az, int N) {
    double *d_Cx, *d_Cy, *d_Cz, *d_mass;
    double *d_ax, *d_ay, *d_az;

    cudaMalloc(&d_Cx, N * sizeof(double));
    cudaMalloc(&d_Cy, N * sizeof(double));
    cudaMalloc(&d_Cz, N * sizeof(double));
    cudaMalloc(&d_mass, N * sizeof(double));
    cudaMalloc(&d_ax, N * sizeof(double));
    cudaMalloc(&d_ay, N * sizeof(double));
    cudaMalloc(&d_az, N * sizeof(double));

    cudaMemcpy(d_Cx, stars->Cx, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Cy, stars->Cy, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Cz, stars->Cz, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mass, stars->mass, N * sizeof(double), cudaMemcpyHostToDevice);

    int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    cuda_kernel<<<numBlocks, BLOCK_SIZE>>>(d_Cx, d_Cy, d_Cz, d_mass, d_ax, d_ay, d_az, N);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
        return;
    }

    cudaMemcpy(ax, d_ax, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(ay, d_ay, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(az, d_az, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_Cx);
    cudaFree(d_Cy);
    cudaFree(d_Cz);
    cudaFree(d_mass);
    cudaFree(d_ax);
    cudaFree(d_ay);
    cudaFree(d_az);
}



