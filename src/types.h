#ifndef TYPES_H
#define TYPES_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define G 4.30091e-6F    // pc^3 / (M_sun * s^2)
#define KAPPA 4.74047F   // Factor de conversión: mas/año a km/s
#define V_GAL 220       //Velocidad media de rotacion galáctica

#define SIGMA 3.0857e13 //Factor de conversion de km/s a parsecs

// Movimiento del Sol respecto al LSR en km/s
#define U_SOL 11.1F   //hacia centro galáctico
#define V_SOL 232.24F //tangencial al eje de rotacion
#define W_SOL 7.25F   //perpendicular a plano galáctico

#define DT 10000          // Paso de tiempo en segundos
#define STEPS 10          // Pasos de la simulacion

typedef struct {
    unsigned long *id;
    float *ra, *dec, *parallax, *pmra, *pmdec, *radial_velocity;
    float *Cx, *Cy, *Cz;
    float *Vx, *Vy, *Vz;
    float *mean_g, *color;
    float *mass;
    int size;
    int capacity;
} Star;
#endif // TYPES_H
