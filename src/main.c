#include "file_handler.h"
#include "types.h"
#include "simulation.h"

// Función para calcular la distancia considerando paralaje negativa o incierta
float calcular_distancia(float paralaje) {
    if (paralaje > 0.1) {
        return 1 / paralaje;
    }
    if (paralaje < 0.1 && paralaje < 0) {
        return 10.0F; // Paralaje muy pequeña, usamos estimación aproximada
    }
    return 4.0F; // Valor bayesiano promedio basado en Bailer-Jones (2018)
}

void complete_data(Star *stars) {
    struct timeval start, end;

    const float R[3][3] = {
        {-0.05487556F, -0.87343709F, -0.48383502F},
        {+0.49410943F, -0.44482963F, +0.74698225F},
        {-0.86766615F, -0.19807637F, +0.45598380F}
    };

    gettimeofday(&start, NULL);

    for (int i = 0; i < stars->size; i++) {
        float d = calcular_distancia(stars->parallax[i]);

        float ra_rad = stars->ra[i] * (M_PI / 180.0);
        float dec_rad = stars->dec[i] * (M_PI / 180.0);

        float cos_ra = cosf(ra_rad), sin_ra = sinf(ra_rad);
        float cos_dec = cosf(dec_rad), sin_dec = sinf(dec_rad);

        float X_E = cos_dec * cos_ra;
        float Y_E = cos_dec * sin_ra;
        float Z_E = sin_dec;

        stars->Cx[i] = (R[0][0] * X_E + R[0][1] * Y_E + R[0][2] * Z_E) * d * 1000;
        stars->Cy[i] = (R[1][0] * X_E + R[1][1] * Y_E + R[1][2] * Z_E) * d * 1000;
        stars->Cz[i] = (R[2][0] * X_E + R[2][1] * Y_E + R[2][2] * Z_E) * d * 1000;

        if (stars->radial_velocity[i] == 0) {
            float C[3] = {stars->Cx[i], stars->Cy[i], stars->Cz[i]};
            float V[3];
            float T[3][2] = {
                {-sin_ra, -cos_ra * sin_dec},
                {cos_ra, -sin_ra * sin_dec},
                {0, cos_dec}
            };

            for (int j = 0; j < 3; j++) {
                V[j] = KAPPA * d * (T[j][0] * stars->pmra[i] + T[j][1] * stars->pmdec[i]);
            }

            float norm_C = sqrtf(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);
            float V_radial_sin = (C[0]*V[0] + C[1]*V[1] + C[2]*V[2]) / norm_C;

            float R_gal = sqrtf(C[0]*C[0] + C[1]*C[1]);
            float V_rot = (V_GAL * C[1]) / R_gal;
            float V_sol_corr = (U_SOL * C[0] + V_SOL * C[1] + W_SOL * C[2]) / norm_C;

            stars->radial_velocity[i] = V_radial_sin + V_rot - V_sol_corr;
        }

        float Vt_ra = KAPPA * stars->pmra[i] * d;
        float Vt_dec = KAPPA * stars->pmdec[i] * d;

        float Vr_eq[3] = {
            stars->radial_velocity[i] * X_E - Vt_ra * sin_ra - Vt_dec * sin_dec * cos_ra,
            stars->radial_velocity[i] * Y_E + Vt_ra * cos_ra - Vt_dec * sin_dec * sin_ra,
            stars->radial_velocity[i] * Z_E + Vt_dec * cos_dec
        };

        stars->Vx[i] = (R[0][0] * Vr_eq[0] + R[0][1] * Vr_eq[1] + R[0][2] * Vr_eq[2]) / SIGMA;
        stars->Vy[i] = (R[1][0] * Vr_eq[0] + R[1][1] * Vr_eq[1] + R[1][2] * Vr_eq[2]) / SIGMA;
        stars->Vz[i] = (R[2][0] * Vr_eq[0] + R[2][1] * Vr_eq[1] + R[2][2] * Vr_eq[2]) / SIGMA;

        float M_G = stars->mean_g[i] * log10f(d) - 5.0F;
        float L_Lsun = powf(10.0F, (4.75F - M_G) / 2.5F);
        float color = stars->color[i];
        float mass;

        if (color < 0.3)
            mass = powf(L_Lsun, 1.0F / 3.8F);
        else if (color < 0.8)
            mass = powf(L_Lsun, 1.0F / 4.0F);
        else if (color < 1.5)
            mass = powf(L_Lsun, 1.0F / 4.5F);
        else
            mass = powf(L_Lsun, 1.0F / 5.0F);

        mass *= 1.0F + 0.1F * (1.0F - expf(-powf(color - 0.8F, 2.0F)));
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
