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
                double dist = sqrt(dist_sq) + 1e-10; // Evitar división por 0
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
    const double epsilon = 1e-10; // Para evitar división por cero (al paralelizar no deberia ser necesario...)
    const __m512d G_vec = _mm512_set1_pd(G); // Almacenamos G
    const __m512d eps_vec = _mm512_set1_pd(epsilon); // Almacenamos Epsilon
    for (int i = 0; i < N; i++) {
        //inicializamos vectores
        __m512d ax_vec = _mm512_setzero_pd();
        __m512d ay_vec = _mm512_setzero_pd();
        __m512d az_vec = _mm512_setzero_pd();
        //cargamos coordenadas de la estrella
        __m512d ix = _mm512_set1_pd(estrellas[i].C[0]);
        __m512d iy = _mm512_set1_pd(estrellas[i].C[1]);
        __m512d iz = _mm512_set1_pd(estrellas[i].C[2]);
        for (int j = 0; j < N; j += 8) { // Procesamos en bloques de 8
            __m512d jx = _mm512_loadu_pd(&estrellas[j].C[0]);
            __m512d jy = _mm512_loadu_pd(&estrellas[j].C[1]);
            __m512d jz = _mm512_loadu_pd(&estrellas[j].C[2]);
            __m512d mass = _mm512_loadu_pd(&estrellas[j].mass);

            __m512d dx = _mm512_sub_pd(jx, ix);
            __m512d dy = _mm512_sub_pd(jy, iy);
            __m512d dz = _mm512_sub_pd(jz, iz);

            __m512d dist_sq = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
            __m512d dist = _mm512_sqrt_pd(_mm512_add_pd(dist_sq, eps_vec));


            __m512d force = _mm512_div_pd(_mm512_mul_pd(G_vec, mass), _mm512_mul_pd(dist_sq, dist));

            ax_vec = _mm512_fmadd_pd(force, dx, ax_vec);
            ay_vec = _mm512_fmadd_pd(force, dy, ay_vec);
            az_vec = _mm512_fmadd_pd(force, dz, az_vec);
        }
        // Devolvemos los resultados a memoria
        double temp_ax[8], temp_ay[8], temp_az[8];
        _mm512_storeu_pd(temp_ax, ax_vec);
        _mm512_storeu_pd(temp_ay, ay_vec);
        _mm512_storeu_pd(temp_az, az_vec);

        for (int k = 0; k < 8 && i + k < N; k++) {
            ax[i + k] += temp_ax[k];
            ay[i + k] += temp_ay[k];
            az[i + k] += temp_az[k];
        }
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
        fprintf(file, "X: %.6f Y = %.6f, Z = %.6f\n", estrellas[i].C[0], estrellas[i].C[1], estrellas[i].C[2]);
    }
}
