#include "file_handler.h"
#include "types.h"
#include "simulation.h"

void print_estrellas(Star *stars) {
    for (int i = 300000; i < 301000 && i < stars->size; i++) {
        printf("------------------------------------------------------------\n");
        printf("ID: %lu\n", stars->id[i]);
        printf("RA: %.4f   DEC: %.4f   Parallax: %.4f   Radial Velocity: %.4f\n",
               stars->ra[i], stars->dec[i], stars->parallax[i], stars->radial_velocity[i]);
        printf("Mean G: %.4f   Color: %.4f   Mass: %.4f\n",
               stars->mean_g[i], stars->color[i], stars->mass[i]);
        printf("Position (X, Y, Z):   (%.20lf, %.20lf, %.20lf)\n",
               stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        printf("Velocity (Vx, Vy, Vz): (%.20lf, %.20lf, %.20lf)\n",
               stars->Vx[i], stars->Vy[i], stars->Vz[i]);
        printf("------------------------------------------------------------\n");
    }
}
int main(int argc, char *argv[]) {
    Star *estrellas = malloc(sizeof(Star));
    memset(estrellas, 0, sizeof(Star));
    if (argc<2) {
        perror("No se ha introducido ningún archivo");
        return -1;
    }
    int num_estrellas = getstarsfromfile(argv[1], estrellas);
    if (num_estrellas < 0) {
        perror("No se encontro ninguna estrella");
        return -1;
    }
    //simulate(estrellas,1000);
    //print_estrellas(estrellas);
    //for (int i = 1000; i <= 1000000; i *= 10) {
    //    simulate(estrellas, i);
    //}
    free_stars(estrellas);
    return 0;
}
