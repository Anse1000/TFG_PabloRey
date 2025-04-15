#include "simulation.h"

//Funcion de cómputo de aceleración
void compute_aceleration(Star *estrellas, double *ax, double *ay, double *az, int N) {
    for (int i = 0; i < N; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                //calcular distancia a la estrella
                double dx = estrellas[j].C[0] - estrellas[i].C[0];
                double dy = estrellas[j].C[1] - estrellas[i].C[1];
                double dz = estrellas[j].C[2] - estrellas[i].C[2];
                double dist_sq = dx * dx + dy * dy + dz * dz;
                double dist = sqrt(dist_sq);
                //calcular fuerza aplicada a la estrella
                double force = -G * estrellas[j].mass / (dist_sq * dist);
                ax[i] += force * dx;
                ay[i] += force * dy;
                az[i] += force * dz;
            }
        }
    }
}

#ifdef AVX_512
void compute_aceleration_avx512(Star *estrellas, double *ax, double *ay, double *az, int N) {
    const __m512d G_vec = _mm512_set1_pd(-G);

    for (int i = 0; i < N; i++) {
        double ax_i = 0.0, ay_i = 0.0, az_i = 0.0;
        // Posición de la estrella i
        __m512d ix = _mm512_set1_pd(estrellas[i].C[0]);
        __m512d iy = _mm512_set1_pd(estrellas[i].C[1]);
        __m512d iz = _mm512_set1_pd(estrellas[i].C[2]);

        for (int j = 0; j < N; j += 8) {
            __mmask8 mask = (j + 7 < N) ? 0xFF : (1 << (N - j)) - 1;

            // Cargar posiciones y masas manualmente
            double temp_jx[8], temp_jy[8], temp_jz[8];
            float temp_mass[8];
            int indices[8];

            for (int k = 0; k < 8 && (j + k) < N; k++) {
                temp_jx[k] = estrellas[j + k].C[0];
                temp_jy[k] = estrellas[j + k].C[1];
                temp_jz[k] = estrellas[j + k].C[2];
                temp_mass[k] = estrellas[j + k].mass;
                indices[k] = j + k;
            }

            __m512d jx = _mm512_maskz_loadu_pd(mask, temp_jx);
            __m512d jy = _mm512_maskz_loadu_pd(mask, temp_jy);
            __m512d jz = _mm512_maskz_loadu_pd(mask, temp_jz);
            __m256 mass_f = _mm256_maskz_loadu_ps(mask, temp_mass);
            __m512d mass = _mm512_cvtps_pd(mass_f);

            // Calcular diferencias
            __m512d dx = _mm512_sub_pd(jx, ix);
            __m512d dy = _mm512_sub_pd(jy, iy);
            __m512d dz = _mm512_sub_pd(jz, iz);

            // Calcular distancia
            __m512d dist_sq = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
            dist_sq = _mm512_add_pd(dist_sq, epsilon);
            __m512d dist = _mm512_sqrt_pd(dist_sq);

            // Crear máscara para i != j
            __m256i index_vec = _mm256_maskz_loadu_epi32(mask, indices);
            __m512i index64 = _mm512_cvtepi32_epi64(index_vec);
            __mmask8 neq_mask = _mm512_cmpneq_epu64_mask(index64, _mm512_set1_epi64(i));

            // Calcular fuerza solo donde i != j
            __m512d force = _mm512_maskz_div_pd(neq_mask,
                                  _mm512_mul_pd(G_vec, mass),
                                  _mm512_mul_pd(dist_sq, dist));

            __m512d fx = _mm512_maskz_mul_pd(neq_mask, force, dx);
            __m512d fy = _mm512_maskz_mul_pd(neq_mask, force, dy);
            __m512d fz = _mm512_maskz_mul_pd(neq_mask, force, dz);

            ax_i += _mm512_reduce_add_pd(fx);
            ay_i += _mm512_reduce_add_pd(fy);
            az_i += _mm512_reduce_add_pd(fz);
        }

        ax[i] += ax_i;
        ay[i] += ay_i;
        az[i] += az_i;
    }
}
#endif

// Función principal de simulación
void simulate(Star *estrellas, const int N) {
    struct timeval start, end;
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
            estrellas[i].V[0] += 0.5 * DT * ax[i];
            estrellas[i].V[1] += 0.5 * DT * ay[i];
            estrellas[i].V[2] += 0.5 * DT * az[i];

            // Actualizar posición
            estrellas[i].C[0] += estrellas[i].V[0] * DT;
            estrellas[i].C[1] += estrellas[i].V[1] * DT;
            estrellas[i].C[2] += estrellas[i].V[2] * DT;
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
            estrellas[i].V[0] += 0.5 * DT * ax[i];
            estrellas[i].V[1] += 0.5 * DT * ay[i];
            estrellas[i].V[2] += 0.5 * DT * az[i];
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
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n",estrellas[i].id ,estrellas[i].C[0], estrellas[i].C[1], estrellas[i].C[2]);
    }
}
