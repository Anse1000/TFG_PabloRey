#include "simulation.h"
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include "aux_fun.h"
#include "octree.h"

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
                 double *ay, double *az, double *seconds)  {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    compute_acceleration_bh(stars, tree, node_idx, index, theta, ax, ay, az);
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

// Función principal de simulación
void simulate(Star *estrellas,const int steps, const long N, const char* outputfile) {
    struct timeval start, end;
    double DT2 = 0.5 * DT;
    double *ax = malloc(N * sizeof(double));
    double *ay = malloc(N * sizeof(double));
    double *az = malloc(N * sizeof(double));

    gettimeofday(&start, NULL);
    printf("Iniciando simulacion usando %d threads\n",omp_get_max_threads()); fflush(stdout);
    for (int step = 0; step < steps; step++) {
        printf("Iniciando paso %d\n", step + 1); fflush(stdout);
        Octree *octree = build_tree(estrellas);
        printf("Iniciando fase 1\n"); fflush(stdout);
        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            compute_acceleration_bh(estrellas, octree, 0, i, 0.2, &ax[i], &ay[i], &az[i]);
            // Leapfrog integration: actualizar velocidad a mitad de paso
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
            // Actualizar posición
            estrellas->Cx[i] = fma(DT, estrellas->Vx[i], estrellas->Cx[i]);
            estrellas->Cy[i] = fma(DT, estrellas->Vy[i], estrellas->Cy[i]);
            estrellas->Cz[i] = fma(DT, estrellas->Vz[i], estrellas->Cz[i]);
        }
        free_tree(octree);
        //Reconstruir con nuevas posiciones
        octree = build_tree(estrellas);
        printf("Iniciando fase 2\n"); fflush(stdout);
        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            compute_acceleration_bh(estrellas, octree, 0, i, 0.2, &ax[i], &ay[i], &az[i]);
            // Completar actualización de velocidad
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
        }
        printf("Paso %d realizado\n", step + 1);
        fflush(stdout);
    }
    gettimeofday(&end, NULL);
    double seconds = get_seconds(start, end);
    int hours = (int)(seconds / 3600);
    int minutes = ((int)seconds % 3600) / 60;
    double remaining_seconds = fmod(seconds, 60.0);
    printf("Simuladas %ld estrellas en %02d:%02d:%05.2f (hh:mm:ss) usando %d threads\n",
       N, hours, minutes, remaining_seconds, omp_get_max_threads());
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
    FILE *file = fopen(outputfile, "w");
    if (!file) {
        fprintf(stderr, "Error al abrir el archivo de salida\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n", estrellas->id[i], estrellas->Cx[i], estrellas->Cy[i],
                estrellas->Cz[i]);
    }
}

void test_simulation(Star *estrellas) {
    double THETA[10]= {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    int indexes[20];
    double seconds[20], seconds_bh[200];
    for (int i = 0; i < 20; i++) {
        indexes[i] = rand() % estrellas->size;
    }
    double ax[20] = {0}, ay[20] = {0}, az[20] = {0};
    double axb[200] = {0}, ayb[200] = {0}, azb[200] = {0};

    Octree *octree = build_tree(estrellas);

    for (int i = 0; i < 20; i++) {
        compute_aceleration_single(estrellas, &ax[i], &ay[i], &az[i], indexes[i], &seconds[i]);
        for (int j = 0; j < 10; j++) {
            int idex=i*10+j;
            aux_time_bh(estrellas, octree, 0, indexes[i], THETA[j], &axb[idex], &ayb[idex], &azb[idex], &seconds_bh[idex]);
        }
        printf("------------------------------------------------------\n");
        printf("Estrella: %d\n", indexes[i]);
        printf("Referencia:              X= %+e Y= %+e Z= %+e  %f segundos\n", ax[i], ay[i], az[i], seconds[i]);
        for (int j = 0; j < 10; j++) {
            int idex=i*10+j;
            printf("BarnesHut THETA %0.1f:     X= %+e Y= %+e Z= %+e  %f segundos\n",THETA[j], axb[idex], ayb[idex], azb[idex],seconds_bh[idex]);
        }
    }

    free_tree(octree);
}
