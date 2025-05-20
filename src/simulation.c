#include "simulation.h"

//Funcion de cómputo de aceleración
void compute_aceleration(Star *estrellas, double *ax, double *ay, double *az,int N) {
    for (int i = 0; i < N; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                //calcular distancia a la estrella
                double dx = estrellas->Cx[i] - estrellas->Cx[j];
                double dy = estrellas->Cy[i] - estrellas->Cy[j];
                double dz = estrellas->Cz[i] - estrellas->Cz[j];
                double dist_sq = dx * dx + dy * dy + dz * dz;
                double dist = sqrt(dist_sq);
                //calcular fuerza aplicada a la estrella
                double force = -G * estrellas->mass[j] / (dist_sq * dist);
                ax[i] = fma(force, dx, ax[i]);
                ay[i] = fma(force, dy, ay[i]);
                az[i] = fma(force, dz, az[i]);
            }
        }
    }
}

#ifdef AVX_512
void compute_aceleration_avx512(const Star *stars, double *ax, double *ay, double *az, int N) {
    const __m512d G_vec = _mm512_set1_pd(-G);

    for (int i = 0; i < N; i++) {
        __m512d xi = _mm512_set1_pd(stars->Cx[i]);
        __m512d yi = _mm512_set1_pd(stars->Cy[i]);
        __m512d zi = _mm512_set1_pd(stars->Cz[i]);

        __m512d sum_ax = _mm512_setzero_pd();
        __m512d sum_ay = _mm512_setzero_pd();
        __m512d sum_az = _mm512_setzero_pd();

        for (int j = 0; j < N; j += 8) {
            // Cargar 8 elementos
            __m512d xj = _mm512_load_pd(&stars->Cx[j]);
            __m512d yj = _mm512_load_pd(&stars->Cy[j]);
            __m512d zj = _mm512_load_pd(&stars->Cz[j]);
            __m512d mj = _mm512_load_pd(&stars->mass[j]);

            __m512d dx = _mm512_sub_pd(xi, xj);
            __m512d dy = _mm512_sub_pd(yi, yj);
            __m512d dz = _mm512_sub_pd(zi, zj);

            __m512d dist2 = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
            __m512d dist = _mm512_sqrt_pd(dist2);
            __m512d dist3 = _mm512_mul_pd(dist2, dist);

            // Fuerza = -G * m / r^3
            __m512d inv = _mm512_div_pd(G_vec, dist3);
            __m512d F = _mm512_mul_pd(inv, mj);

            // Crear máscara para i != j
            __mmask8 mask = _mm512_cmpneq_epi64_mask(
                _mm512_set1_epi64(i),
                _mm512_set_epi64(j+7, j+6, j+5, j+4, j+3, j+2, j+1, j+0)
            );

            // sum += F * dX con máscara
            sum_ax = _mm512_mask3_fmadd_pd(F, dx, sum_ax, mask);
            sum_ay = _mm512_mask3_fmadd_pd(F, dy, sum_ay, mask);
            sum_az = _mm512_mask3_fmadd_pd(F, dz, sum_az, mask);
        }

        // Reducción horizontal para obtener el valor final
        ax[i] = _mm512_reduce_add_pd(sum_ax);
        ay[i] = _mm512_reduce_add_pd(sum_ay);
        az[i] = _mm512_reduce_add_pd(sum_az);
    }
}
#endif

// Función principal de simulación
void simulate(Star *estrellas, const int N) {
    struct timeval start, end;
    double DT2 = 0.5 * DT;
    double *ax = malloc(N * sizeof(double));
    double *ay = malloc(N * sizeof(double));
    double *az = malloc(N * sizeof(double));
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
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
            // Actualizar posición
            estrellas->Cx[i] = fma(DT, estrellas->Vx[i], estrellas->Cx[i]);
            estrellas->Cy[i] = fma(DT, estrellas->Vy[i], estrellas->Cy[i]);
            estrellas->Cz[i] = fma(DT,estrellas->Vz[i], estrellas->Cz[i]);
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
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
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
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n",estrellas->id[i] ,estrellas->Cx[i], estrellas->Cy[i], estrellas->Cz[i]);
    }
}
