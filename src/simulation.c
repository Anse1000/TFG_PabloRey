#include "simulation.h"
//Funcion de cómputo de aceleración
void compute_aceleration(Star *estrellas, float *ax, float *ay, float *az, int N) {
    for (int i = 0; i < N; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                //calcular distancia a la estrella
                float dx = estrellas->Cx[i] - estrellas->Cx[j];
                float dy = estrellas->Cy[i] - estrellas->Cy[j];
                float dz = estrellas->Cz[i] - estrellas->Cz[j];
                float dist_sq = dx * dx + dy * dy + dz * dz;
                float dist = sqrtf(dist_sq);
                //calcular fuerza aplicada a la estrella
                float force = -G * estrellas->mass[j] / (dist_sq * dist);
                ax[i] += force * dx;
                ay[i] += force * dy;
                az[i] += force * dz;
            }
        }
    }
}

#ifdef AVX_512
void compute_aceleration_avx512(const Star *stars, float *ax, float *ay, float *az, int N) {
    const __m512 G_vec = _mm512_set1_ps(-G);

    for (int i = 0; i < N; i++) {
        __m512 xi = _mm512_set1_ps(stars->Cx[i]);
        __m512 yi = _mm512_set1_ps(stars->Cy[i]);
        __m512 zi = _mm512_set1_ps(stars->Cz[i]);

        __m512 sum_ax = _mm512_setzero_ps();
        __m512 sum_ay = _mm512_setzero_ps();
        __m512 sum_az = _mm512_setzero_ps();

        for (int j = 0; j < N; j += 16) {
            // Cargar 8 elementos
            __m512 xj = _mm512_load_ps(&stars->Cx[j]);
            __m512 yj = _mm512_load_ps(&stars->Cy[j]);
            __m512 zj = _mm512_load_ps(&stars->Cz[j]);
            __m512 mj = _mm512_load_ps(&stars->mass[j]);

            __m512 dx = _mm512_sub_ps(xi, xj);
            __m512 dy = _mm512_sub_ps(yi, yj);
            __m512 dz = _mm512_sub_ps(zi, zj);

            __m512 dist2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
            __m512 dist = _mm512_sqrt_ps(dist2);
            __m512 dist3 = _mm512_mul_ps(dist2, dist);

            // Fuerza = -G * m / r^3
            __m512 inv = _mm512_div_ps(G_vec, dist3);
            __m512 F = _mm512_mul_ps(inv, mj);

            __mmask16 mask = _mm512_cmpneq_epi32_mask(
                _mm512_set1_epi32(i),
                _mm512_set_epi32(
                    j + 15, j + 14, j + 13, j + 12, j + 11, j + 10, j + 9, j + 8,
                    j + 7, j + 6, j + 5, j + 4, j + 3, j + 2, j + 1, j + 0
                )
            );

            // sum += F * dX con máscara
            sum_ax = _mm512_mask3_fmadd_ps(F, dx, sum_ax, mask);
            sum_ay = _mm512_mask3_fmadd_ps(F, dy, sum_ay, mask);
            sum_az = _mm512_mask3_fmadd_ps(F, dz, sum_az, mask);
        }

        // Reducción horizontal para obtener el valor final
        ax[i] = _mm512_reduce_add_ps(sum_ax);
        ay[i] = _mm512_reduce_add_ps(sum_ay);
        az[i] = _mm512_reduce_add_ps(sum_az);
    }
}
#endif

// Función principal de simulación
void simulate(Star *estrellas, const int N) {
    struct timeval start, end;
    float *ax = malloc(N * sizeof(float));
    float *ay = malloc(N * sizeof(float));
    float *az = malloc(N * sizeof(float));
    gettimeofday(&start, NULL);
    for (int step = 0; step < STEPS; step++) {
#ifdef AVX_512
        compute_aceleration_avx512(estrellas, ax, ay, az, N);
#elif CUDA
        compute_aceleration_CUDA(estrellas, ax, ay, az, N);
#else
        compute_aceleration(estrellas, ax, ay, az, N);

#endif
        for (int i = 0; i < N; i++) {
            // Leapfrog integration: actualizar velocidad a mitad de paso
            estrellas->Vx[i] += 0.5F * DT * ax[i];
            estrellas->Vy[i] += 0.5F * DT * ay[i];
            estrellas->Vz[i] += 0.5F * DT * az[i];

            // Actualizar posición
            estrellas->Cx[i] += estrellas->Vx[i] * DT;
            estrellas->Cy[i] += estrellas->Vy[i] * DT;
            estrellas->Cz[i] += estrellas->Vz[i] * DT;
        }
#ifdef AVX_512
        compute_aceleration_avx512(estrellas, ax, ay, az, N);
#elif CUDA
        compute_aceleration_CUDA(estrellas, ax, ay, az, N);
#else
        compute_aceleration(estrellas, ax, ay, az, N);
#endif
        for (int i = 0; i < N; i++) {
            // Completar actualización de velocidad
            estrellas->Vx[i] += 0.5F * DT * ax[i];
            estrellas->Vy[i] += 0.5F * DT * ay[i];
            estrellas->Vz[i] += 0.5F * DT * az[i];
        }
        printf("Paso %d realizado\n", step + 1);
        fflush(stdout);
    }
    gettimeofday(&end, NULL);
    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Simuladas %d estrellas en %.2f segundos\n", N, seconds);
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
    FILE *file = fopen("galactic_positions.txt", "w");
    if (!file) {
        fprintf(stderr, "Error al abrir el archivo\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n", estrellas->id[i], estrellas->Cx[i], estrellas->Cy[i],
                estrellas->Cz[i]);
    }
}
