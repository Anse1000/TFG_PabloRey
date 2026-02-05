#include "simulation.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "../aux_fun.h"
#include "octree_cpu.h"
#include "../file_handler.h"

// --- Halo NFW
double halo_accel(double r, double *ax, double *ay, double *az,
                  double dx, double dy, double dz) {
    // concentración c = R200/rs ~ 10, R200 ~ 200 kpc
    double x = r / rs;
    double f = log(1.0 + x) - x / (1.0 + x);
    double Menc = M200 * f / (log(1.0 + 10.0) - 10.0 / 11.0);

    double acc = -G * Menc / (r * r * r);
    *ax = fma(acc, dx, *ax);
    *ay = fma(acc, dy, *ay);
    *az = fma(acc, dz, *az);
    return acc;
}

// --- Bulbo Hernquist
double bulge_accel(double r, double *ax, double *ay, double *az,
                   double dx, double dy, double dz) {
    double acc = -G * MBULGE / ((r + A) * (r + A)) / r;
    *ax = fma(acc, dx, *ax);
    *ay = fma(acc, dy, *ay);
    *az = fma(acc, dz, *az);
    return acc;
}

// --- Wrapper total (halo + bulbo)
void analytic_accel(double x, double y, double z,
                    double *ax, double *ay, double *az) {
    double dx = x;
    double dy = y;
    double dz = z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    if (r > 0) {
        halo_accel(r, ax, ay, az, dx, dy, dz);
        bulge_accel(r, ax, ay, az, dx, dy, dz);
    }
}

//funcion de prueba: calcula la aceleracion de una sola estrella con TODAS
void compute_aceleration_single(const Star *stars, double *ax, double *ay, double *az, const unsigned long index,
                                double *seconds) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    for (unsigned long i = 0; i < stars->size; i++) {
        if (i != index) {
            //calcular distancia a la estrella
            double dx = stars->Cx[i] - stars->Cx[index];
            double dy = stars->Cy[i] - stars->Cy[index];
            double dz = stars->Cz[i] - stars->Cz[index];
            double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;
            double dist = sqrt(dist_sq);
            //calcular fuerza aplicada a la estrella
            double force = -G * stars->mass[i] / (dist_sq * dist);
            *ax = fma(force, dx, *ax);
            *ay = fma(force, dy, *ay);
            *az = fma(force, dz, *az);
        }
    }
    analytic_accel(stars->Cx[index], stars->Cy[index], stars->Cz[index], ax, ay, az);
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

void compute_acceleration_bh(const Star *stars, const Octree *tree,
                             long node_idx, long star_idx, double theta,
                             double *ax, double *ay, double *az) {
    // Ignorar si es una hoja con la misma estrella
    if (tree->star_index[node_idx] == star_idx)
        return;

    // Diferencia de posición
    double dx = tree->com_x[node_idx] - stars->Cx[star_idx];
    double dy = tree->com_y[node_idx] - stars->Cy[star_idx];
    double dz = tree->com_z[node_idx] - stars->Cz[star_idx];
    double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;
    double dist = sqrt(dist_sq);

    double s = 2.0 * tree->half_size[node_idx]; // ancho total del nodo

    if ((s / dist) < theta || tree->star_index[node_idx] >= 0) {
        // Tratar como nodo lejano o hoja
        double force = -G * tree->mass[node_idx] / (dist_sq * dist);
        *ax = fma(force, dx, *ax);
        *ay = fma(force, dy, *ay);
        *az = fma(force, dz, *az);
    } else {
        // Recursión en hijos
        for (int i = 0; i < 8; i++) {
            long child = tree->children[node_idx][i];
            if (child != INVALID_INDEX) {
                compute_acceleration_bh(stars, tree, child, star_idx, theta, ax, ay, az);
            }
        }
    }
}

void aux_time_bh(const Star *stars, const Octree *tree, long node_idx, long index, double theta, double *ax,
                 double *ay, double *az, double *seconds) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    compute_acceleration_bh(stars, tree, node_idx, index, theta, ax, ay, az);
    analytic_accel(stars->Cx[index], stars->Cy[index], stars->Cz[index], ax, ay, az);
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

// Función principal de simulación
void simulate(Star *estrellas, const int steps, const long N, const char *outputfile, float DT) {
    struct timeval start, end;
    double *ax = malloc(N * sizeof(double));
    double *ay = malloc(N * sizeof(double));
    double *az = malloc(N * sizeof(double));
    float cx, cy, cz;
    float hs, min_node_size;
    gettimeofday(&start, NULL);
    printf("******************************************************\n");
    printf("Iniciando simulacion con %d threads\n", omp_get_max_threads());
    fflush(stdout);
    float DT2 = 0.5F * DT;
    write_results_hdf5(estrellas, outputfile, "cpu_results", -1);
    printf("Simulando %d pasos de %.0f años (Total: %.0f años)\n", steps, DT * 1000000, DT * steps * 1000000);
    printf("******************************************************\n");
    fflush(stdout);
    // Inicializar árbol y aceleraciones antes del bucle principal
    compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
    Octree *octree = build_tree(estrellas, cx, cy, cz, hs, min_node_size);

    // Array de aceleraciones (reutilizable entre pasos)
#pragma omp parallel for
    for (long i = 0; i < N; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
        compute_acceleration_bh(estrellas, octree, 0, i, THETA, &ax[i], &ay[i], &az[i]);
        analytic_accel(estrellas->Cx[i], estrellas->Cy[i], estrellas->Cz[i], &ax[i], &ay[i], &az[i]);
    }

    for (int step = 0; step < steps; step++) {
        struct timeval step_start, step_end;
        printf("****************** Iniciando paso %d ******************\n", step + 1);
        gettimeofday(&step_start, NULL);

        // === HalfKick usando aceleraciones previas ===
        printf("\t  Fase 1: HalfKick-Drift\n");
        fflush(stdout);
#pragma omp parallel for
        for (long i = 0; i < N; i++) {
            // half-kick
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);

            // drift
            estrellas->Cx[i] = fma(DT, estrellas->Vx[i], estrellas->Cx[i]);
            estrellas->Cy[i] = fma(DT, estrellas->Vy[i], estrellas->Cy[i]);
            estrellas->Cz[i] = fma(DT, estrellas->Vz[i], estrellas->Cz[i]);
        }

        // Liberar y reconstruir árbol con nuevas posiciones
        free_tree(octree);
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
        octree = build_tree(estrellas, cx, cy, cz, hs, min_node_size);

        // === Calcular nuevas aceleraciones ===
        printf("\t  Fase 2: HalfKick\n");
        fflush(stdout);
#pragma omp parallel for
        for (long i = 0; i < N; i++) {
            ax[i] = ay[i] = az[i] = 0.0;
            compute_acceleration_bh(estrellas, octree, 0, i, THETA, &ax[i], &ay[i], &az[i]);
            analytic_accel(estrellas->Cx[i], estrellas->Cy[i], estrellas->Cz[i], &ax[i], &ay[i], &az[i]);

            // completar el segundo half-kick con las nuevas aceleraciones
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
        }

        // Guardar resultados
        write_results_hdf5(estrellas, outputfile, "cpu_results", step);

        gettimeofday(&step_end, NULL);
        double step_seconds = get_seconds(step_start, step_end);
        printf("************ Paso %d finalizado en %6.0f segundos ************\n",
               step + 1, step_seconds);
        fflush(stdout);
    }

    gettimeofday(&end, NULL);
    double seconds = get_seconds(start, end);
    int hours = (int) (seconds / 3600);
    int minutes = ((int) seconds % 3600) / 60;
    double remaining_seconds = fmod(seconds, 60.0);
    printf("Simulacion de %ld estrellas completada en %02d:%02d:%05.2f (hh:mm:ss)\n",
           N, hours, minutes, remaining_seconds);
    printf("Resultados guardados en %s\n", outputfile);
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
}

void test_simulation(Star *estrellas) {
    double THETAS[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    int indexes[20];
    double seconds[20], seconds_bh[200];
    for (int i = 0; i < 20; i++) {
        indexes[i] = rand() % estrellas->size;
    }
    double ax[20] = {0}, ay[20] = {0}, az[20] = {0};
    double axb[200] = {0}, ayb[200] = {0}, azb[200] = {0};
    float cx, cy, cz;
    float hs, min_node_size;
    compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size,MIN_SUBDIVISIONS);
    Octree *octree = build_tree(estrellas, cx, cy, cz, hs, min_node_size);

    for (int i = 0; i < 20; i++) {
        compute_aceleration_single(estrellas, &ax[i], &ay[i], &az[i], indexes[i], &seconds[i]);
        for (int j = 0; j < 10; j++) {
            int idex = i * 10 + j;
            aux_time_bh(estrellas, octree, 0, indexes[i], THETAS[j], &axb[idex], &ayb[idex], &azb[idex],
                        &seconds_bh[idex]);
        }
        printf("------------------------------------------------------\n");
        printf("Estrella: %d\n", indexes[i]);
        printf("Position (X, Y, Z):   (%.20lf, %.20lf, %.20lf)\n",
               estrellas->Cx[i], estrellas->Cy[i], estrellas->Cz[i]);
        printf("Velocity (Vx, Vy, Vz): (%.20lf, %.20lf, %.20lf)\n",
               estrellas->Vx[i], estrellas->Vy[i], estrellas->Vz[i]);
        printf("Referencia:              X= %.20lf Y= %.20lf Z= %.20lf  %f segundos\n", ax[i], ay[i], az[i],
               seconds[i]);
        for (int j = 0; j < 10; j++) {
            int idex = i * 10 + j;
            printf("BarnesHut THETA %0.1f:     X= %.20lf Y= %.20lf Z= %.20lf  %f segundos\n", THETAS[j], axb[idex],
                   ayb[idex],
                   azb[idex], seconds_bh[idex]);
        }
    }
    free_tree(octree);
}
// --- FUNCIONES DE APOYO PARA TESTS ---

// Crea una copia reducida de las estrellas para tests rápidos
Star* clone_subsample(const Star *original, int factor) {
    long N_sub = original->size / factor;
    Star *sub = (Star*)malloc(sizeof(Star));
    sub->size = N_sub;
    sub->id = (unsigned long*)malloc(N_sub * sizeof(long));
    sub->Cx = (double*)malloc(N_sub * sizeof(double)); sub->Cy = (double*)malloc(N_sub * sizeof(double)); sub->Cz = (double*)malloc(N_sub * sizeof(double));
    sub->Vx = (double*)malloc(N_sub * sizeof(double)); sub->Vy = (double*)malloc(N_sub * sizeof(double)); sub->Vz = (double*)malloc(N_sub * sizeof(double));
    sub->mass = (float*)malloc(N_sub * sizeof(float));

    for(long i=0; i < N_sub; i++) {
        long idx = i * factor;
        sub->id[i] = original->id[idx];
        sub->Cx[i] = original->Cx[idx]; sub->Cy[i] = original->Cy[idx]; sub->Cz[i] = original->Cz[idx];
        sub->Vx[i] = original->Vx[idx]; sub->Vy[i] = original->Vy[idx]; sub->Vz[i] = original->Vz[idx];
        sub->mass[i] = original->mass[idx] * (float)factor; // Escalado de masa
    }
    printf("Muestra reducida\n");
    return sub;
}

void free_temp_stars(Star *s) {
    free(s->id); free(s->Cx); free(s->Cy); free(s->Cz);
    free(s->Vx); free(s->Vy); free(s->Vz); free(s->mass);
    free(s);
}

// --- LOS 3 TESTS DE CPU ---

// 1. TEST UNITARIO: Órbita Circular (Validación física de constantes y aceleración analítica)
void run_cpu_unit_orbit_test() {
    printf("\n[TEST 1] Órbita Circular Unitaria (Física y Constantes)\n");
    Star *s = (Star*)malloc(sizeof(Star));
    memset(s, 0, sizeof(Star));
    s->size = 1;
    s->Cx = (double*)malloc(sizeof(double)); s->Cy = (double*)malloc(sizeof(double)); s->Cz = (double*)malloc(sizeof(double));
    s->Vx = (double*)malloc(sizeof(double)); s->Vy = (double*)malloc(sizeof(double)); s->Vz = (double*)malloc(sizeof(double));
    s->mass = (float*)malloc(sizeof(float));

    // 1. Posición inicial
    s->Cx[0] = 8.2; s->Cy[0] = 0.0; s->Cz[0] = 0.0;
    s->mass[0] = 1.0f;

    // 2. Calcular aceleración en el punto inicial para obtener la v_circular exacta
    double ax0=0, ay0=0, az0=0;
    analytic_accel(s->Cx[0], s->Cy[0], s->Cz[0], &ax0, &ay0, &az0);
    double a_mag = sqrt(ax0*ax0 + ay0*ay0 + az0*az0);

    // v = sqrt(r * a) para órbita circular
    double v_mag = sqrt(8.2 * a_mag);
    s->Vx[0] = 0.0; s->Vy[0] = v_mag; s->Vz[0] = 0.0;

    printf("  Aceleración inicial: %.6le kpc/Myr^2\n", a_mag);
    printf("  Velocidad circular calculada: %.6f kpc/Myr (aprox %.2f km/s)\n", v_mag, v_mag / 0.0010227);

    double r_init = s->Cx[0];
    const float dt = 0.1f;
    const float dt2 = dt * 0.5f;

    for(int i=0; i<1000; i++) {
        double ax=0, ay=0, az=0;
        analytic_accel(s->Cx[0], s->Cy[0], s->Cz[0], &ax, &ay, &az);
        s->Vx[0] = fma(dt2, ax, s->Vx[0]); s->Vy[0] = fma(dt2, ay, s->Vy[0]); s->Vz[0] = fma(dt2, az, s->Vz[0]);
        s->Cx[0] = fma(dt, s->Vx[0], s->Cx[0]); s->Cy[0] = fma(dt, s->Vy[0], s->Cy[0]); s->Cz[0] = fma(dt, s->Vz[0], s->Cz[0]);
        ax=ay=az=0;
        analytic_accel(s->Cx[0], s->Cy[0], s->Cz[0], &ax, &ay, &az);
        s->Vx[0] = fma(dt2, ax, s->Vx[0]); s->Vy[0] = fma(dt2, ay, s->Vy[0]); s->Vz[0] = fma(dt2, az, s->Vz[0]);
    }

    double r_final = sqrt(s->Cx[0]*s->Cx[0] + s->Cy[0]*s->Cy[0] + s->Cz[0]*s->Cz[0]);
    printf("  Radio tras 100 Myr: Inicial=%.4f kpc, Final=%.4f kpc\n", r_init, r_final);
    printf("  Estabilidad: %s\n", (fabs(r_init-r_final) < 0.1) ? "PASADO" : "FALLADO (revisar G o unidades)");
    free_temp_stars(s);
}

// 2. TEST ESCALABILIDAD: Analítico Masivo (Stress-test de OpenMP y RAM con 1.100M)
void run_cpu_massive_analytical_test(const Star *original) {
    printf("\n[TEST 2] Analítico Masivo (N=%lu - Stress-test OpenMP)\n", original->size);
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Creamos vectores de aceleración temporales para no tocar los del main
    double *tax = (double*)malloc(original->size * sizeof(double));
    double *tay = (double*)malloc(original->size * sizeof(double));
    double *taz = (double*)malloc(original->size * sizeof(double));

    printf("  Calculando aceleraciones analíticas para todo el dataset...\n");
    #pragma omp parallel for
    for (size_t i = 0; i < original->size; i++) {
        tax[i] = tay[i] = taz[i] = 0.0;
        analytic_accel(original->Cx[i], original->Cy[i], original->Cz[i], &tax[i], &tay[i], &taz[i]);
    }

    gettimeofday(&end, NULL);
    printf("  Completado en %.2f segundos (%.2f millones estrellas/seg)\n",
            get_seconds(start, end), (original->size/1e6)/get_seconds(start, end));
    free(tax); free(tay); free(taz);
}
// Calcula el momento angular total en el eje Z: Lz = sum( m * (x*vy - y*vx) )
double calculate_total_Lz(const Star *s) {
    double total_Lz = 0.0;
#pragma omp parallel for reduction(+:total_Lz)
    for (size_t i = 0; i < s->size; i++) {
        total_Lz += (double)s->mass[i] *
                    (s->Cx[i] * s->Vy[i] - s->Cy[i] * s->Vx[i]);
    }
    return total_Lz;
}

// Helper: Ejecuta un paso completo de integración KDK (idéntico al de la simulación real)
void perform_kdk_step(Star *s, Octree **tree, float DT) {
    float DT2 = DT * 0.5f;
    long N = s->size;
    float cx, cy, cz, hs, min_node_size;

    // 1. Primer Half-Kick + Drift
    #pragma omp parallel for
    for(long i=0; i<N; i++) {
        double ax=0, ay=0, az=0;
        compute_acceleration_bh(s, *tree, 0, i, THETA, &ax, &ay, &az);
        analytic_accel(s->Cx[i], s->Cy[i], s->Cz[i], &ax, &ay, &az);

        s->Vx[i] += ax * DT2; s->Vy[i] += ay * DT2; s->Vz[i] += az * DT2;
        s->Cx[i] += s->Vx[i] * DT; s->Cy[i] += s->Vy[i] * DT; s->Cz[i] += s->Vz[i] * DT;
    }

    // 2. Reconstrucción del árbol
    free_tree(*tree);
    compute_root_bounds(s, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
    *tree = build_tree(s, cx, cy, cz, hs, min_node_size);

    // 3. Segundo Half-Kick con las nuevas aceleraciones
    #pragma omp parallel for
    for(long i=0; i<N; i++) {
        double ax=0, ay=0, az=0;
        compute_acceleration_bh(s, *tree, 0, i, THETA, &ax, &ay, &az);
        analytic_accel(s->Cx[i], s->Cy[i], s->Cz[i], &ax, &ay, &az);
        s->Vx[i] += ax * DT2; s->Vy[i] += ay * DT2; s->Vz[i] += az * DT2;
    }
}

// 3. TEST ALGORÍTMICO: Reversibilidad con esquema real
void run_cpu_reversibility_subsample(const Star *original, float DT_sim) {
    int factor = 100;
    printf("\n[TEST 3] Reversibilidad Barnes-Hut (Full KDK, DT=%.2f, Muestra 1/%d)\n", DT_sim, factor);

    Star *sub = clone_subsample(original, factor);
    long N = sub->size;
    double *orig_x = (double*)malloc(N * sizeof(double));
    double *orig_y = (double*)malloc(N * sizeof(double));
    double *orig_z = (double*)malloc(N * sizeof(double));
    for(long i=0; i<N; i++) { orig_x[i] = sub->Cx[i]; orig_y[i] = sub->Cy[i]; orig_z[i] = sub->Cz[i]; }

    float cx, cy, cz, hs, min_node_size;
    compute_root_bounds(sub, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
    Octree *tree = build_tree(sub, cx, cy, cz, hs, min_node_size);

    printf("  Evolucionando 2 pasos ADELANTE...\n");
    for(int s=0; s<2; s++) perform_kdk_step(sub, &tree, DT_sim);

    printf("  Evolucionando 2 pasos ATRÁS (DT = %.2f)...\n", -DT_sim);
    for(int s=0; s<2; s++) perform_kdk_step(sub, &tree, -DT_sim);

    double err = 0;
    #pragma omp parallel for reduction(+:err)
    for(long i=0; i<N; i++) {
        err += sqrt(pow(sub->Cx[i]-orig_x[i],2) + pow(sub->Cy[i]-orig_y[i],2) + pow(sub->Cz[i]-orig_z[i],2));
    }

    printf("  Error final MAE: %.6le kpc\n", err/N);
    printf("  Resultado: %s\n", (err/N < 1e-4) ? "PASADO" : "REVISAR PRECISIÓN");

    free_tree(tree);
    free(orig_x); free(orig_y); free(orig_z);
    free_temp_stars(sub);
}

// 4. TEST DE SISTEMAS: Conservación del Momento Angular (Lz)
void run_cpu_angular_momentum_test(const Star *original, float DT_sim) {
    int factor = 100;
    printf("\n[TEST 4] Conservación del Momento Angular Lz (DT=%.2f, Muestra 1/%d)\n", DT_sim, factor);

    Star *sub = clone_subsample(original, factor);
    double Lz_init = calculate_total_Lz(sub);

    float cx, cy, cz, hs, min_node_size;
    compute_root_bounds(sub, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
    Octree *tree = build_tree(sub, cx, cy, cz, hs, min_node_size);

    printf("  Evolucionando 5 pasos para medir deriva...\n");
    for(int s=0; s<5; s++) perform_kdk_step(sub, &tree, DT_sim);

    double Lz_final = calculate_total_Lz(sub);
    double rel_err = fabs(Lz_final - Lz_init) / (fabs(Lz_init) + 1e-20);

    printf("  Variación relativa Lz: %.4le\n", rel_err);
    free_tree(tree);
    free_temp_stars(sub);
}

void run_full_cpu_validation(Star *estrellas) {
    float DT = 1.0f;
    printf("\n#######################################################\n");
    printf("      INICIANDO SUITE DE VALIDACIÓN COMPLETA (CPU)\n");
    printf("#######################################################\n");

    run_cpu_unit_orbit_test();
    run_cpu_reversibility_subsample(estrellas, DT);
    run_cpu_massive_analytical_test(estrellas);
    run_cpu_angular_momentum_test(estrellas, DT);

    printf("\n#######################################################\n");
}
void test_theta(Star *stars, int num_samples) {
    double thetas[] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                       0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    int num_thetas = 20;

    // Arrays de resultados
    double sum_errors[20] = {0}, sum_sq_errors[20] = {0};
    double max_errors[20] = {0}, sum_bh_times[20] = {0};
    double total_ref_time = 0;

    float cx, cy, cz, hs, min_node_size;
    compute_root_bounds(stars, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
    Octree *tree = build_tree(stars, cx, cy, cz, hs, min_node_size);
    printf("Iniciando %d muestras con %d hilos...\n", num_samples, omp_get_max_threads()); fflush(stdout);

    #pragma omp parallel
    {
        // Semilla local por hilo para evitar contención en rand()
        unsigned int seed = time(NULL) ^ omp_get_thread_num();

        #pragma omp for reduction(+:total_ref_time, sum_errors[:20], sum_sq_errors[:20], sum_bh_times[:20]) \
                       reduction(max:max_errors[:20])
        for (int k = 0; k < num_samples; k++) {
            unsigned long idx = rand_r(&seed) % stars->size;

            // 1. Referencia (Lento: ~8s)
            double ref_ax = 0, ref_ay = 0, ref_az = 0, ref_time = 0;
            compute_aceleration_single(stars, &ref_ax, &ref_ay, &ref_az, idx, &ref_time);
            double ref_mag = sqrt(ref_ax*ref_ax + ref_ay*ref_ay + ref_az*ref_az);
            total_ref_time += ref_time;

            // 2. Test de Barnes-Hut (Rápido: <0.04s)
            for (int t = 0; t < num_thetas; t++) {
                double bh_ax = 0, bh_ay = 0, bh_az = 0, bh_time = 0;
                aux_time_bh(stars, tree, 0, idx, thetas[t], &bh_ax, &bh_ay, &bh_az, &bh_time);

                double dx = bh_ax - ref_ax, dy = bh_ay - ref_ay, dz = bh_az - ref_az;
                double rel_error = (ref_mag > 0.0) ? (sqrt(dx*dx + dy*dy + dz*dz) / ref_mag) : 0.0;

                sum_errors[t] += rel_error;
                sum_sq_errors[t] += (rel_error * rel_error);
                sum_bh_times[t] += bh_time;
                if (rel_error > max_errors[t]) max_errors[t] = rel_error;
            }
        }
    }

    // Informe final (promediado sobre num_samples)
    printf("\n--- RESULTADOS FINALES (%d muestras) ---\n", num_samples);
    double reference_time = total_ref_time / num_samples;
    printf("Tiempo medio Fuerza bruta: %.6f\n",reference_time);
    for (int t = 0; t < num_thetas; t++) {
        double m_err = sum_errors[t] / num_samples;
        double std_dev = sqrt(fmax(0, (sum_sq_errors[t] / num_samples) - (m_err * m_err)));
        printf("Theta: %.2f | Err Med: %.2e | Err Max: %.2e | StdDev: %.2e | T_medio: %.6f\n",
                thetas[t], m_err, max_errors[t], std_dev, sum_bh_times[t] / num_samples);
    }
    free_tree(tree);
}
void benchmark_tree_construction(Star *estrellas, int iteraciones) {
    if (iteraciones < 2) iteraciones = 2; // Mínimo 2 para tener warm-up

    double *tiempos = (double*) malloc(iteraciones * sizeof(double));

    printf("\n============================================================\n");
    printf(" BENCHMARK: Construcción del Árbol Barnes-Hut\n");
    printf(" Estrellas: %lu | Hilos OpenMP: %d | Iteraciones: %d\n",
           estrellas->size, omp_get_max_threads(), iteraciones);
    printf("============================================================\n");
    fflush(stdout);

    // Variables para bounds (se recalculan en cada iteración para ser realistas)
    float cx, cy, cz, hs, min_node_size;

    for (int i = 0; i < iteraciones; i++) {

        // Sincronización previa para que todos los hilos arranquen a la vez
        #pragma omp barrier
        double start = omp_get_wtime();

        // --- INICIO FASE CRÍTICA ---
        compute_root_bounds(estrellas, &cx, &cy, &cz, &hs, &min_node_size, MIN_SUBDIVISIONS);
        Octree *tree = build_tree(estrellas, cx, cy, cz, hs, min_node_size);
        // ---------------------------

        double end = omp_get_wtime();
        tiempos[i] = end - start;

        printf("   Iteración %02d: %.6f s %s\n",
               i, tiempos[i], (i == 0) ? "(Warm-up - Descartada)" : "");

        free_tree(tree);
    }

    // --- CÁLCULO ESTADÍSTICO ---
    double suma = 0.0, suma_sq = 0.0;
    int conteo_valido = iteraciones - 1;

    // Empezamos desde i=1 para ignorar la primera pasada (cache fría)
    for (int i = 1; i < iteraciones; i++) {
        suma += tiempos[i];
    }
    double media = suma / conteo_valido;

    for (int i = 1; i < iteraciones; i++) {
        suma_sq += pow(tiempos[i] - media, 2);
    }
    double desviacion = sqrt(suma_sq / conteo_valido);

    printf("------------------------------------------------------------\n");
    printf(" RESULTADO FINAL (Media +/- Desviación):\n");
    printf(" Tiempo: %.6f s +/- %.6f s\n", media, desviacion);
    printf(" Rate:   %.2f Millones de estrellas/seg\n", (estrellas->size / 1e6) / media);
    printf("============================================================\n");

    free(tiempos);
}