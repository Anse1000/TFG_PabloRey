#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/time.h>

#define G 4.30091e-6  // pc^3 / (M_sun * s^2)
#define KAPPA 4.74047   // Factor de conversión: mas/año a km/s
#define V_GAL 220       //Velocidad media de rotacion galáctica
#define DT 100000          // Paso de tiempo en segundos
#define SIGMA 3.0857e13  //Factor de conversion de km/s a parsecs

// Movimiento del Sol respecto al LSR en km/s
#define U_SOL 11.1 //hacia centro galáctico
#define V_SOL 232.24 //tangencial al eje de rotacion
#define W_SOL 7.25 //perpendicular a plano galáctico

#define STEPS 10

#define MAX_LINE 2000
#define DELIMITER ","

#define TYPE_UL  0  // unsigned long
#define TYPE_D   1  // double
#define TYPE_F   2  // float

typedef struct {
    unsigned long id;
    double ra, dec, parallax, pmra, pmdec, radial_velocity, mass;
    double C[3], V[3];
    float mean_g, color;
} Star;

typedef struct {
    void *ptr;
    int type;
} FieldMap;

int read_file(char *filename, Star **estrellas, int *N, int *allocated_size) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        return -1;
    }

    char line[MAX_LINE];

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == 'i') continue;

        char *token = strtok(line, DELIMITER);
        int is_valid = 1;
        char *tokens[21];
        int i = 0;

        while (token && i < 21) {
            if ((i == 3 || i == 7 || i == 8) && token[0] == 'n') {
                is_valid = 0;
                break;
            }
            tokens[i++] = token;
            token = strtok(NULL, DELIMITER);
        }

        if (!is_valid) continue;

        //Se reserva de 10000 en 10000 para evitar tantas llamadas a realloc
        if (*N >= *allocated_size) {
            *allocated_size += 10000;  // Aumentar de a 10000 elementos
            Star *temp = realloc(*estrellas, *allocated_size * sizeof(Star));
            if (!temp) {
                fclose(file);
                return -1;
            }
            *estrellas = temp;
        }

        Star *current = &(*estrellas)[*N];
        (*N)++;

        const FieldMap fields[] = {
                {&current->id,                    TYPE_UL},
                {&current->ra,                    TYPE_D},
                {&current->dec,                   TYPE_D},
                {&current->parallax,              TYPE_D},
                {&current->pmra,                  TYPE_D},
                {&current->pmdec,                 TYPE_D},
                {&current->radial_velocity,       TYPE_D},
                {&current->mean_g,                TYPE_F},
                {&current->color,                 TYPE_F}
        };

        for (int j = 0; j < 9; j++) {
            switch (fields[j].type) {
                case TYPE_UL:
                    *(unsigned long *) fields[j].ptr = strtoul(tokens[j], NULL, 10);
                    break;
                case TYPE_D:
                    *(double *) fields[j].ptr = strtod(tokens[j], NULL);
                    break;
                case TYPE_F:
                    *(float *) fields[j].ptr = strtof(tokens[j], NULL);
                    break;
            }
        }
    }
    fclose(file);
    return 0;
}

int getstarsfromfile(char *dirname,Star **estrellas) {
    struct timeval start,end;
    char path[1000];
    struct dirent *filedir;
    int dirname_length=strlen(dirname);
    int N=0,count=0,allocated=0;
    int status=0;
    DIR *dir = opendir(dirname);
    if (dir == NULL) {
        perror("No se pudo abrir el directorio");
        return -1;
    }
    strcpy(path, dirname);
    strcat(path, "/");

    gettimeofday(&start, NULL);
    while ((filedir = readdir(dir)) != NULL) {
        if (filedir->d_name[0] == '.') continue;
        strcpy(path + dirname_length + 1, filedir->d_name);
        status=read_file(path,estrellas,&N,&allocated);
        if (status!=0) {
            printf("Error en la lectura de archivos\n");
            return -1;
        }
        count++;
        if (count%10==0) {
            printf("Leidos %d archivos...\n",count);
            fflush(stdout);
        }

    }
    //Ajustar memoria a tamaño exacto
    Star *temp = realloc(*estrellas, N * sizeof(Star));  // Redimensionar al número exacto de estrellas leídas
    if (temp != NULL) {
        *estrellas = temp;
    } else {
        perror("No se pudo redimensionar la memoria.\n");
    }
    gettimeofday(&end, NULL);
    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("\nLeídas y trasladadas %d estrellas a memoria ocupando %.2f MB en %.2f segundos\n", N,
       N * sizeof(Star) / (1024.0 * 1024.0),seconds);
    fflush(stdout);
    fflush(0);
    return N;
}

//Funcion de cómputo de aceleración
void compute_acceleration(Star *estrellas, double *ax, double *ay, double *az, int N) {
    for (int i = 0; i < N; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                //calcular distancia a la estrella
                double dx = estrellas[j].C[0] - estrellas[i].C[0];
                double dy = estrellas[j].C[1] - estrellas[i].C[1];
                double dz = estrellas[j].C[2] - estrellas[i].C[2];
                double dist_sq = dx * dx + dy * dy + dz * dz;
                double dist = sqrt(dist_sq) + 1e-10; // Evitar división por 0
                //calcular fuerza aplicada a la estrella
                double force = G * estrellas[j].mass / (dist_sq * dist);
                ax[i] += force * dx;
                ay[i] += force * dy;
                az[i] += force * dz;
            }
        }
    }
}

// Función principal de simulación
void simulate(Star *estrellas, const int N) {
    struct timeval start,end;
    double *ax = malloc(N*sizeof(double));
    double *ay = malloc(N*sizeof(double));
    double *az = malloc(N*sizeof(double));
    gettimeofday(&start, NULL);
    for (int step = 0; step < STEPS; step++) {
        compute_acceleration(estrellas, ax, ay, az, N);
        for (int i = 0; i < N; i++) {
            // Leapfrog integration: actualizar velocidad a mitad de paso
            estrellas[i].V[0] += 0.5 * DT * ax[i];
            estrellas[i].V[1] += 0.5 * DT * ay[i];
            estrellas[i].V[2] += 0.5 * DT * az[i];

            // Actualizar posición
            estrellas[i].C[0] += estrellas[i].V[0] * DT;
            estrellas[i].C[1] += estrellas[i].V[1] * DT;
            estrellas[i].C[2] += estrellas[i].V[2] * DT;
        }
        compute_acceleration(estrellas, ax, ay, az, N);
        for (int i = 0; i < N; i++) {
            // Completar actualización de velocidad
            estrellas[i].V[0] += 0.5 * DT * ax[i];
            estrellas[i].V[1] += 0.5 * DT * ay[i];
            estrellas[i].V[2] += 0.5 * DT * az[i];
        }
        printf("Paso %d realizado\n", step+1);
        fflush(stdout);
    }
    gettimeofday(&end, NULL);
    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Simuladas %d estrellas en %.2f segundos\n",N,seconds);
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
    FILE *file = fopen("galactic_positions.txt", "w");
    if (!file) {
        fprintf(stderr, "Error al abrir el archivo\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "X: %.6f Y = %.6f, Z = %.6f\n", estrellas[i].C[0], estrellas[i].C[1], estrellas[i].C[2]);
    }
}

// Función para calcular la distancia considerando paralaje negativa o incierta
double calcular_distancia(double paralaje) {
    if (paralaje > 0.1) {
        return 1 / paralaje;
    }
    if (paralaje<0.1 && paralaje<0) {
        return 10.0;  // Paralaje muy pequeña, usamos estimación aproximada
    }
    return 4.0;  // Valor bayesiano promedio basado en Bailer-Jones (2018)
}

void complete_data(Star *estrellas, int N) {
    struct timeval start,end;
    // Matriz de transformación ecuatorial a galáctico
    const double R[3][3] = {
        {-0.05487556, -0.87343709, -0.48383502},
        {+0.49410943, -0.44482963, +0.74698225},
        {-0.86766615, -0.19807637, +0.45598380}
    };
    gettimeofday(&start, NULL);
    for (int i = 0; i < N; i++) {
        double d = calcular_distancia(estrellas[i].parallax);
        // Convertir grados a radianes
        double ra_rad = estrellas[i].ra * (M_PI / 180.0);
        double dec_rad = estrellas[i].dec * (M_PI / 180.0);
        //Precálculos para mejorar eficiencia
        double cos_ra = cos(ra_rad);
        double sin_ra = sin(ra_rad);
        double cos_dec = cos(dec_rad);
        double sin_dec = sin(dec_rad);
        /***Paso 1: Coordenadas galácticas X,Y,Z almacenadas en estrella->C ***/
        //Convertir a coordenadas cartesianas ecuatoriales
        double X_E = cos_dec * cos_ra;
        double Y_E = cos_dec * sin_ra;
        double Z_E = sin_dec;
        // Transformación a coordenadas galácticas cartesianas
        for (int j = 0; j < 3; ++j) {
            estrellas[i].C[j] = (R[j][0] * X_E + R[j][1] * Y_E + R[j][2] * Z_E)*d*1000;
        }
        /***Si no se incluye la velocidad radial se realiza una estimación***/
        if (estrellas[i].radial_velocity == 0) {
            double C[3] = {estrellas[i].C[0], estrellas[i].C[1],
                           estrellas[i].C[2]}; //variable para facilitar escritura
            double V[3];
            // Matriz de transformación de proper motion en RA DEC a vx, vy, vz
            double T[3][2] = {
                    {-sin_ra, -cos_ra * sin_dec},
                    {cos_ra,  -sin_ra * sin_dec},
                    {0,            cos_dec}
            };
            // Multiplicación de matriz por vector
            for (int j = 0; j < 3; j++) {
                V[j] = KAPPA * d * (T[j][0] * estrellas[i].pmra + T[j][1] * estrellas[i].pmdec);
            }
            // Estimación de la velocidad radial
            double V_radial_sin =
                    (C[0] * V[0] + C[1] * V[1] + C[2] * V[2]) / sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
            // Distancia radial desde el centro galáctico
            double R_gal = sqrt(C[0] * C[0] + C[1] * C[1]);
            // Corregir por la rotación galáctica
            double V_rot = (V_GAL * C[1]) / R_gal;
            // Corregir por la velocidad del Sol
            double V_sol_correccion =
                    (U_SOL * C[0] + V_SOL * C[1] + W_SOL * C[2]) / sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);

            // Resultado final con corrección por rotación y velocidad del Sol
            estrellas[i].radial_velocity = V_radial_sin + V_rot - V_sol_correccion;
        }
        /*** Paso 2: conseguir vectores de velocidad completos Vx,Vy,Vz ***/
        // Calcular velocidad tangencial
        double Vt_ra = KAPPA * estrellas[i].pmra * d;
        double Vt_dec = KAPPA * estrellas[i].pmdec * d;

        // Vector de velocidad en el sistema ecuatorial
        double Vr_eq[3] = {
                estrellas[i].radial_velocity * X_E - Vt_ra * sin_ra - Vt_dec * sin_dec * cos_ra,
                estrellas[i].radial_velocity * Y_E + Vt_ra * cos_ra - Vt_dec * sin_dec * sin_ra,
                estrellas[i].radial_velocity * Z_E + Vt_dec * cos_dec
        };

        // Multiplicar la matriz de transformación por el vector de velocidad
        for (int j = 0; j < 3; ++j) {
            estrellas[i].V[j] = (R[j][0] * Vr_eq[0] + R[j][1] * Vr_eq[1] + R[j][2] * Vr_eq[2])/SIGMA;
        }
        /*** Paso 3: calcular masa ***/
        // Calcular magnitud absoluta M_G
        double M_G = estrellas[i].mean_g * log10(d) - 5.0;

        // Calcular luminosidad relativa L/L_sun
        double L_Lsun = pow(10.0, (4.75 - M_G) / 2.5);

        // Estimar la masa en función de la luminosidad y el color
        double mass;
        if (estrellas[i].color < 0.3) {
            // Estrellas calientes (tipo O, B, A)
            mass = pow(L_Lsun, 1.0 / 3.8);
        } else if (estrellas[i].color < 0.8) {
            // Estrellas tipo F-G (similares al Sol)
            mass = pow(L_Lsun, 1.0 / 4.0);
        } else if (estrellas[i].color < 1.5) {
            // Estrellas tipo K (más frías y menos masivas)
            mass = pow(L_Lsun, 1.0 / 4.5);
        } else {
            // Enanas rojas (M)
            mass = pow(L_Lsun, 1.0 / 5.0);
        }
        // Corrección adicional basada en el color para mejorar precisión
        mass *= 1.0 + 0.1 * (1.0 - exp(-pow(estrellas[i].color - 0.8, 2.0)));

        //Masa basada en masas solares
        estrellas[i].mass = mass;
    }
    gettimeofday(&end, NULL);
    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Completados datos de %d estrellas en %lf segundos\n",N, seconds);
}

void print_estrellas(Star *estrellas) {
    for (int i = 300000; i < 301000; i++) {//Se escogen estrellas del medio para no ver siempre las mismas
        printf("------------------------------------------------------------\n");
        printf("ID: %lu\n", estrellas[i].id);
        printf("RA: %.4f   DEC: %.4f   Parallax: %.4f   Radial Velocity: %.4f\n",
               estrellas[i].ra, estrellas[i].dec, estrellas[i].parallax, estrellas[i].radial_velocity);
        printf("Mean G: %.4f   Color: %.4f   Mass: %.4f\n",
               estrellas[i].mean_g, estrellas[i].color, estrellas[i].mass);
        printf("Position (X, Y, Z):   (%.4f, %.4f, %.4f)\n",
               estrellas[i].C[0], estrellas[i].C[1], estrellas[i].C[2]);
        printf("Velocity (Vx, Vy, Vz): (%.20f, %.20f, %.20f)\n",
               estrellas[i].V[0], estrellas[i].V[1], estrellas[i].V[2]);
        printf("------------------------------------------------------------\n");
    }
}

int main() {
    Star *estrellas = NULL;
    int num_estrellas = getstarsfromfile("reducidos",&estrellas);
    if (num_estrellas < 0) {
        perror("No se encontro ninguna estrella");
    }
    complete_data(estrellas,num_estrellas);
    //print_estrellas(estrellas);
    for (int i = 1000; i <= 100000; i*=10) {
        simulate(estrellas,i);
    }
    free(estrellas);  // Liberar memoria
    return 0;
}
