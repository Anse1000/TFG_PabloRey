#ifndef TYPES_H
#define TYPES_H
#include <stdio.h>

#define G 4.498540135e-12  // Constante gravitatoria universal
#define KAPPA 4.74047   // AS/año * parsecs -> km/s
#define V_GAL 220       // Velocidad media de rotacion galáctica
#define SIGMA 0.0010227 // Factor de conversion de km/s a kiloparsecs/Milones de años

// Movimiento del Sol respecto al LSR en km/s
#define U_SOL 11.1   // Hacia centro galáctico
#define V_SOL 12.24  // Tangencial al eje de rotacion
#define W_SOL 7.25   // Perpendicular a plano galáctico

#define EPSILON 0.00001 //Valor para evitar divisiones por cero o muy cercanas a cero
#define MIN_SUBDIVISIONS 1e7 //Numero maximo de subdivisiones del arbol
#define THETA 0.2 //Apertura del arbol BH

// parámetros de halo
#define M200 1.0e12   // masa virial en Msolares
#define rs 20.0     // radio de escala en kpc
// parámetros bulbo
#define MBULGE 1.0e10  // Msolares
#define A 0.7     // kpc (escala bulbo)

typedef struct {
    unsigned long *id;
    double *ra, *dec, *distance, *pmra, *pmdec, *radial_velocity;
    double *Cx, *Cy, *Cz;
    double *Vx, *Vy, *Vz;
    float *color;
    float *mass;
    size_t size;
    size_t capacity;
} Star;


#endif // TYPES_H
