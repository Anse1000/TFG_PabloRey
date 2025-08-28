#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "aux_fun.h"
#include "file_handler.h"
#include "types.h"
#include "simulation.h"
#ifdef CUDA
#include "cuda_functions.cuh"
#endif

void print_estrellas(Star *stars) {
    for (size_t i = 300000; i < 301000 && i < stars->size; i++) {
        printf("------------------------------------------------------------\n");
        printf("ID: %lu\n", stars->id[i]);
        printf("RA: %.4f   DEC: %.4f   Distance: %.4f   Radial Velocity: %.4f\n",
               stars->ra[i], stars->dec[i], stars->distance[i], stars->radial_velocity[i]);
        printf("Mean G: %.4f   Color: %.4f   Mass: %.4f\n",
               stars->mean_g[i], stars->color[i], stars->mass[i]);
        printf("Position (X, Y, Z):   (%.20lf, %.20lf, %.20lf)\n",
               stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        printf("Velocity (Vx, Vy, Vz): (%.20lf, %.20lf, %.20lf)\n",
               stars->Vx[i], stars->Vy[i], stars->Vz[i]);
        printf("------------------------------------------------------------\n");
    }
}

void findminmax(Star *stars) {
    double min_rv = stars->radial_velocity[0], max_rv = stars->radial_velocity[0];
    double min_mass = stars->mass[0], max_mass = stars->mass[0];
    double min_dist = stars->distance[0], max_dist = stars->distance[0];
    double min_cx = stars->Cx[0], max_cx = stars->Cx[0];
    double min_cy = stars->Cy[0], max_cy = stars->Cy[0];
    double min_cz = stars->Cz[0], max_cz = stars->Cz[0];

    for (size_t i = 0; i < stars->size; i++) {
        if (stars->radial_velocity[i] < min_rv) min_rv = stars->radial_velocity[i];
        else if (stars->radial_velocity[i] > max_rv) max_rv = stars->radial_velocity[i];

        if (stars->mass[i] < min_mass) min_mass = stars->mass[i];
        else if (stars->mass[i] > max_mass) max_mass = stars->mass[i];

        if (stars->distance[i] < min_dist) min_dist = stars->distance[i];
        else if (stars->distance[i] > max_dist) max_dist = stars->distance[i];

        if (stars->Cx[i] < min_cx) min_cx = stars->Cx[i];
        else if (stars->Cx[i] > max_cx) max_cx = stars->Cx[i];
        if (stars->Cy[i] < min_cy) min_cy = stars->Cy[i];
        else if (stars->Cy[i] > max_cy) max_cy = stars->Cy[i];
        if (stars->Cz[i] < min_cz) min_cz = stars->Cz[i];
        else if (stars->Cz[i] > max_cz) max_cz = stars->Cz[i];
    }

    printf("Rangos de valores:\n");
    printf("Velocidad Radial: %.8f - %.8f\n", min_rv, max_rv);
    printf("Masa: %.8f - %.8f\n", min_mass, max_mass);
    printf("Distancia: %.8f - %.8f\n", min_dist, max_dist);
    printf("Posición X: %.8f - %.8f\n", min_cx, max_cx);
    printf("Posición Y: %.8f - %.8f\n", min_cy, max_cy);
    printf("Posición Z: %.8f - %.8f\n", min_cz, max_cz);
    fflush(stdout);
}

void write_initial_positions(Star *estrellas, const char *directory) {
    float cx,cy,cz;
    float hs,min_node_size;
    unsigned int offsets[8];
    compute_root_bounds(estrellas,&cx,&cy,&cz,&hs,&min_node_size,MIN_SUBDIVISIONS);
    reorder_stars(estrellas,cx,cy,cz,offsets);
    write_chunks(estrellas,"initial_positions",directory,25,1);
    printf("Posiciones iniciales guardadas en %s\n",directory);
    fflush(stdout);
}

int main(int argc, char *argv[]) {
    int steps = 0;
    Star *estrellas = malloc(sizeof(Star));
    memset(estrellas, 0, sizeof(Star));
    if (argc < 4) {
        fprintf(stderr, "Error: Número incorrecto de argumentos\n");
        fprintf(stderr, "Uso: %s <archivo_estrellas> <archivo_posiciones_iniciales> <archivo_salida> STEPS\n", argv[0]);
        fprintf(stderr, "  <archivo_estrellas>: Carpeta con los datos de entrada de las estrellas\n");
        fprintf(stderr, "  <archivo_posiciones_iniciales>: Carpeta con las posiciones iniciales de las estrellas\n");
        fprintf(stderr, "  <archivo_salida>: Carpeta donde se guardarán los resultados\n");
        fprintf(stderr, "  Numero de pasos de la simulacion\n");
        free(estrellas);
        return -1;
    }
    unsigned long num_estrellas = getstarsfromfile(argv[1], estrellas);
    if (num_estrellas <= 0) {
        perror("No se encontro ninguna estrella");
        return -1;
    }
     write_initial_positions(estrellas,argv[2]);
     free_aux(estrellas);
    steps = atoi(argv[4]);
    if (steps <= 0) {
        perror("Numero de pasos incorrecto");
        return -1;
    }
#ifdef CUDA
    simulate_multi_gpu_unified(estrellas,steps,estrellas->size,argv[3]);
#else
    simulate(estrellas,steps,estrellas->size,argv[3]);
#endif
    //test_simulation(estrellas);
    //test_tree(estrellas);
    free_stars(estrellas);
    return 0;
}
