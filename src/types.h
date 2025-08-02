#ifndef TYPES_H
#define TYPES_H
#include <stdio.h>

#define G 4.30091e-6
#define KAPPA 4.74047   // AS/año * parsecs -> km/s
#define V_GAL 220       // Velocidad media de rotacion galáctica
#define SIGMA 1.0227 // Factor de conversion de km/s a kiloparsecs/Milones de años
#define LOGG_SOL 4.437F  // Logaritmo de la gravedad del SOL

// Movimiento del Sol respecto al LSR en km/s
#define U_SOL 11.1   // Hacia centro galáctico
#define V_SOL 12.24  // Tangencial al eje de rotacion
#define W_SOL 7.25   // Perpendicular a plano galáctico

#define DT 0.0001    // Paso de tiempo en millones de años
#define EPSILON 0.000001
#define MIN_SUBDIVISIONS 1e8

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


#endif // TYPES_H
