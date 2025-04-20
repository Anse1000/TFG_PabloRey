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
                ax[i] = fma(force, dx, ax[i]);
                ay[i] = fma(force, dy, ay[i]);
                az[i] = fma(force, dz, az[i]);
            }
        }
    }
}


extern "C" void compute_aceleration_CUDA(Star *stars, double *ax, double *ay, double *az, int N) {
    double *d_Cx, *d_Cy, *d_Cz, *d_mass;
    double *d_ax, *d_ay, *d_az;
	const size_t size = N * sizeof(double);


	// Usar streams para solapar transferencias
    cudaStream_t stream1, stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    cudaMalloc(&d_Cx, size);
    cudaMalloc(&d_Cy, size);
    cudaMalloc(&d_Cz, size);
    cudaMalloc(&d_mass, size);
    cudaMalloc(&d_ax, size);
    cudaMalloc(&d_ay, size);
    cudaMalloc(&d_az, size);

	// Transferir datos usando streams
    cudaMemcpyAsync(d_Cx, stars->Cx, size, cudaMemcpyHostToDevice, stream1);
    cudaMemcpyAsync(d_Cy, stars->Cy, size, cudaMemcpyHostToDevice, stream1);
    cudaMemcpyAsync(d_Cz, stars->Cz, size, cudaMemcpyHostToDevice, stream2);
    cudaMemcpyAsync(d_mass, stars->mass, size, cudaMemcpyHostToDevice, stream2);

    int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    cuda_kernel<<<numBlocks, BLOCK_SIZE>>>(d_Cx, d_Cy, d_Cz, d_mass, d_ax, d_ay, d_az, N);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
        return;
    }

	// Transferir resultados de vuelta
    cudaMemcpyAsync(ax, d_ax, size, cudaMemcpyDeviceToHost, stream1);
    cudaMemcpyAsync(ay, d_ay, size, cudaMemcpyDeviceToHost, stream2);
    cudaMemcpyAsync(az, d_az, size, cudaMemcpyDeviceToHost, stream2);

    // Sincronizar y limpiar
    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);


    cudaFree(d_Cx);
    cudaFree(d_Cy);
    cudaFree(d_Cz);
    cudaFree(d_mass);
    cudaFree(d_ax);
    cudaFree(d_ay);
    cudaFree(d_az);
}



