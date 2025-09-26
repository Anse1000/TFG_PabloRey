#include "simulation.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "../aux_fun.h"
#include "octree_cpu.h"

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
        write_results(estrellas, outputfile, "cpu_results", step);

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
