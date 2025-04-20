#include "file_handler.h"
#include "types.h"
#include "simulation.h"

// Función para calcular la distancia considerando paralaje negativa o incierta
double calcular_distancia(double paralaje) {
    if (paralaje > 0.1) {
        return 1 / paralaje;
    }
    if (paralaje < 0.1 && paralaje < 0) {
        return 10.0; // Paralaje muy pequeña, usamos estimación aproximada
    }
    return 4.0; // Valor bayesiano promedio basado en Bailer-Jones (2018)
}

void complete_data(Star *stars) {
    struct timeval start, end;

    const double R[3][3] = {
        {-0.05487556, -0.87343709, -0.48383502},
        {+0.49410943, -0.44482963, +0.74698225},
        {-0.86766615, -0.19807637, +0.45598380}
    };

    gettimeofday(&start, NULL);

    for (int i = 0; i < stars->size; i++) {
        double d = calcular_distancia(stars->parallax[i]);

        double ra_rad = stars->ra[i] * (M_PI / 180.0);
        double dec_rad = stars->dec[i] * (M_PI / 180.0);

        double cos_ra = cos(ra_rad), sin_ra = sin(ra_rad);
        double cos_dec = cos(dec_rad), sin_dec = sin(dec_rad);

        double X_E = cos_dec * cos_ra;
        double Y_E = cos_dec * sin_ra;
        double Z_E = sin_dec;

        stars->Cx[i] = (R[0][0] * X_E + R[0][1] * Y_E + R[0][2] * Z_E) * d * 1000;
        stars->Cy[i] = (R[1][0] * X_E + R[1][1] * Y_E + R[1][2] * Z_E) * d * 1000;
        stars->Cz[i] = (R[2][0] * X_E + R[2][1] * Y_E + R[2][2] * Z_E) * d * 1000;

        if (stars->radial_velocity[i] == 0) {
            double C[3] = {stars->Cx[i], stars->Cy[i], stars->Cz[i]};
            double V[3];
            double T[3][2] = {
                {-sin_ra, -cos_ra * sin_dec},
                {cos_ra, -sin_ra * sin_dec},
                {0, cos_dec}
            };

            for (int j = 0; j < 3; j++) {
                V[j] = KAPPA * d * (T[j][0] * stars->pmra[i] + T[j][1] * stars->pmdec[i]);
            }

            double norm_C = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);
            double V_radial_sin = (C[0]*V[0] + C[1]*V[1] + C[2]*V[2]) / norm_C;

            double R_gal = sqrt(C[0]*C[0] + C[1]*C[1]);
            double V_rot = (V_GAL * C[1]) / R_gal;
            double V_sol_corr = (U_SOL * C[0] + V_SOL * C[1] + W_SOL * C[2]) / norm_C;

            stars->radial_velocity[i] = V_radial_sin + V_rot - V_sol_corr;
        }

        double Vt_ra = KAPPA * stars->pmra[i] * d;
        double Vt_dec = KAPPA * stars->pmdec[i] * d;

        double Vr_eq[3] = {
            stars->radial_velocity[i] * X_E - Vt_ra * sin_ra - Vt_dec * sin_dec * cos_ra,
            stars->radial_velocity[i] * Y_E + Vt_ra * cos_ra - Vt_dec * sin_dec * sin_ra,
            stars->radial_velocity[i] * Z_E + Vt_dec * cos_dec
        };

        stars->Vx[i] = (R[0][0] * Vr_eq[0] + R[0][1] * Vr_eq[1] + R[0][2] * Vr_eq[2]) / SIGMA;
        stars->Vy[i] = (R[1][0] * Vr_eq[0] + R[1][1] * Vr_eq[1] + R[1][2] * Vr_eq[2]) / SIGMA;
        stars->Vz[i] = (R[2][0] * Vr_eq[0] + R[2][1] * Vr_eq[1] + R[2][2] * Vr_eq[2]) / SIGMA;

        double M_G = stars->mean_g[i] * log10(d) - 5.0;
        double L_Lsun = pow(10.0, (4.75 - M_G) / 2.5);
        double color = stars->color[i];
        double mass;

        if (color < 0.3)
            mass = pow(L_Lsun, 1.0 / 3.8);
        else if (color < 0.8)
            mass = pow(L_Lsun, 1.0 / 4.0);
        else if (color < 1.5)
            mass = pow(L_Lsun, 1.0 / 4.5);
        else
            mass = pow(L_Lsun, 1.0 / 5.0);

        mass *= 1.0 + 0.1 * (1.0 - exp(-pow(color - 0.8, 2.0)));
        stars->mass[i] = mass;
    }

    gettimeofday(&end, NULL);
    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    printf("Completados datos de %d estrellas en %lf segundos\n", stars->size, seconds);
}

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
    complete_data(estrellas);
    //simulate(estrellas,1000);
    //print_estrellas(estrellas);
    for (int i = 1000; i <= 1000000; i *= 10) {
        simulate(estrellas, i);
    }
    free(estrellas); // Liberar memoria
    return 0;
}
