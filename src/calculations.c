#include "calculations.h"

//Matriz de cambio de coordenadas celestes a galácticas X,Y,Z
//Usando Norte galáctico Epoch 2016.0
const double R[3][3] = {
    {-0.05487556, -0.87343709, -0.48383502},
    {+0.49410943, -0.44482963, +0.74698225},
    {-0.86766615, -0.19807637, +0.45598380}
};

void calculate_mass(double *mass, float gravity, float radius, float mean_g, float color, double d) {
    if (gravity != 0 || radius != 0) {
        // log10(M / M_sun) = logg - logg_sun + 2 * log10(R / R_sun)
        double log_mass = gravity - LOGG_SOL + 2.0 * log10(radius);
        *mass = pow(10.0, log_mass);
    } else {
        double M_G = mean_g - 5.0 * log10(d) + 5.0; //magnitud absoluta a 10 parsecs
        // 2. Relación empírica: log(M) ≈ a * (BP-RP) + b * M_G + c
        // Ajustada para secuencia principal solar-metalicidad (aproximado)
        const double a = -0.15;
        const double b = -0.10;
        const double c = 1.2;

        double log_mass = a * color + b * M_G + c;
        *mass = pow(10.0, log_mass);
    }
}

void calculate_radialvelocity(double *radial_velocity, double sin_ra, double cos_ra, double sin_dec, double cos_dec,
                              double d, double Cx, double Cy, double Cz, double pmra, double pmdec) {
    double V[3];
    double T[3][2] = {
        {-sin_ra, -cos_ra * sin_dec},
        {cos_ra, -sin_ra * sin_dec},
        {0, cos_dec}
    };

    for (int j = 0; j < 3; j++) {
        V[j] = KAPPA * d * (T[j][0] * pmra + T[j][1] * pmdec);
    }

    double norm_C = sqrt(Cx * Cx + Cy * Cy + Cz * Cz);
    double V_radial_sin = (Cx * V[0] + Cy * V[1] + Cz * V[2]) / norm_C;

    double R_gal = sqrt(Cx * Cx + Cy * Cy);
    double V_rot = (V_GAL * Cy) / R_gal;
    double V_sol_corr = (U_SOL * Cx + V_SOL * Cy + W_SOL * Cz) / norm_C;

    *radial_velocity = V_radial_sin + V_rot - V_sol_corr;
}

void calculate_coords(double *Cx, double *Cy, double *Cz, double X_E, double Y_E, double Z_E, double d) {
    *Cx = (R[0][0] * X_E + R[0][1] * Y_E + R[0][2] * Z_E) * d;
    *Cy = (R[1][0] * X_E + R[1][1] * Y_E + R[1][2] * Z_E) * d;
    *Cz = (R[2][0] * X_E + R[2][1] * Y_E + R[2][2] * Z_E) * d;
}

void calculate_vectors(double *Vx, double *Vy, double *Vz, double X_E, double Y_E, double Z_E, double pmra,
                       double pmdec, double rv, double d, double sin_ra, double sin_dec, double cos_ra,
                       double cos_dec) {
    double Vt_ra = KAPPA * pmra * d;
    double Vt_dec = KAPPA * pmdec * d;

    double Vr_eq[3] = {
        rv * X_E - Vt_ra * sin_ra - Vt_dec * sin_dec * cos_ra,
        rv * Y_E + Vt_ra * cos_ra - Vt_dec * sin_dec * sin_ra,
        rv * Z_E + Vt_dec * cos_dec
    };

    *Vx = (R[0][0] * Vr_eq[0] + R[0][1] * Vr_eq[1] + R[0][2] * Vr_eq[2]) / SIGMA;
    *Vy = (R[1][0] * Vr_eq[0] + R[1][1] * Vr_eq[1] + R[1][2] * Vr_eq[2]) / SIGMA;
    *Vz = (R[2][0] * Vr_eq[0] + R[2][1] * Vr_eq[1] + R[2][2] * Vr_eq[2]) / SIGMA;
}

void complete_data(Star *stars) {
    for (int i = 0; i < stars->size; i++) {
        double ra_rad = stars->ra[i] * (M_PI / 180.0);
        double dec_rad = stars->dec[i] * (M_PI / 180.0);

        double cos_ra = cos(ra_rad), sin_ra = sin(ra_rad);
        double cos_dec = cos(dec_rad), sin_dec = sin(dec_rad);

        double X_E = cos_dec * cos_ra;
        double Y_E = cos_dec * sin_ra;
        double Z_E = sin_dec;

        calculate_coords(&stars->Cx[i], &stars->Cy[i], &stars->Cz[i], X_E, Y_E, Z_E, stars->distance[i]);

        if (stars->radial_velocity[i] == 0) {
            calculate_radialvelocity(stars->radial_velocity, sin_ra, cos_ra, sin_dec, cos_dec, stars->distance[i],
                                     stars->Cx[i], stars->Cy[i], stars->Cz[i], stars->pmra[i], stars->pmdec[i]);
        }
        if (stars->mass[i] == 0) {
            calculate_mass(&stars->mass[i], stars->gravity[i], stars->radius[i], stars->mean_g[i], stars->color[i],
                           stars->distance[i]);
        }
        calculate_vectors(stars->Vx, stars->Vy, stars->Vz, X_E, Y_E, Z_E, stars->pmra[i], stars->pmdec[i],
                          stars->radial_velocity[i], stars->distance[i], sin_ra, sin_dec, cos_ra, cos_dec);
    }
}
