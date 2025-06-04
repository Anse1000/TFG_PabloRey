#ifndef TYPES_H
#define TYPES_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define G 4.30091e-6
#define KAPPA 4.74047   // AS/año * parsecs -> km/s
#define V_GAL 220       // Velocidad media de rotacion galáctica
#define SIGMA 1.0227 // Factor de conversion de km/s a kiloparsecs/Milones de años
#define LOGG_SOL 4.437F  // Logaritmo de la gravedad del SOL

// Movimiento del Sol respecto al LSR en km/s
#define U_SOL 11.1   // Hacia centro galáctico
#define V_SOL 12.24  // Tangencial al eje de rotacion
#define W_SOL 7.25   // Perpendicular a plano galáctico

#define DT 10000     // Paso de tiempo en segundos
#define STEPS 10     // Pasos de la simulacion
#define EPSILON 0.000001

typedef struct {
    unsigned long *id;
    double *ra, *dec, *distance, *pmra, *pmdec, *radial_velocity;
    double *Cx, *Cy, *Cz;
    double *Vx, *Vy, *Vz;
    float *mean_g, *color,*radius,*gravity;
    float *mass;
    size_t size;
    size_t capacity;
} Star;

typedef struct {
    double *cx, *cy, *cz;      // centro del cubo
    float *half_size;       // mitad del tamaño del cubo

    float *mass;            // masa total dentro del nodo
    float *com_x, *com_y, *com_z; // centro de masa

    long *star_index;           // -1 si no tiene estrella
    long *first_child_index;  // -1 si es hoja/hijos contiguos: [i, i+8]
    size_t capacity;
    size_t size;
} Octree;

#endif // TYPES_H
