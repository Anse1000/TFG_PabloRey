#include "file_handler.h"

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

static void *safe_realloc(void *ptr, size_t size) {
    void *tmp = realloc(ptr, size);
    if (!tmp) {
        perror("Fallo al hacer el resize de memoria");
        exit(EXIT_FAILURE);
    }
    return tmp;
}

void free_stars(Star *stars) {
    free(stars->id);
    free(stars->ra);
    free(stars->dec);
    free(stars->parallax);
    free(stars->pmdec);
    free(stars->pmra);
    free(stars->radial_velocity);
    free(stars->mean_g);
    free(stars->color);
    free(stars->Cx);
    free(stars->Cy);
    free(stars->Cz);
    free(stars->Vx);
    free(stars->Vy);
    free(stars->Vz);
    free(stars->mass);
    free(stars);
}

void realloc_stars(Star *stars) {
    stars->id = safe_realloc(stars->id, sizeof(long) * stars->capacity);
    stars->ra = safe_realloc(stars->ra, sizeof(double) * stars->capacity);
    stars->dec = safe_realloc(stars->dec, sizeof(double) * stars->capacity);
    stars->parallax = safe_realloc(stars->parallax, sizeof(double) * stars->capacity);
    stars->pmdec = safe_realloc(stars->pmdec, sizeof(double) * stars->capacity);
    stars->pmra = safe_realloc(stars->pmra, sizeof(double) * stars->capacity);
    stars->pmdec = safe_realloc(stars->pmdec, sizeof(double) * stars->capacity);
    stars->radial_velocity = safe_realloc(stars->radial_velocity, sizeof(double) * stars->capacity);
    stars->mean_g = safe_realloc(stars->mean_g, sizeof(float) * stars->capacity);
    stars->color = safe_realloc(stars->color, sizeof(float) * stars->capacity);
    stars->Cx = safe_realloc(stars->Cx, sizeof(double) * stars->capacity);
    stars->Cy = safe_realloc(stars->Cy, sizeof(double) * stars->capacity);
    stars->Cz = safe_realloc(stars->Cz, sizeof(double) * stars->capacity);
    stars->Vx = safe_realloc(stars->Vx, sizeof(double) * stars->capacity);
    stars->Vy = safe_realloc(stars->Vy, sizeof(double) * stars->capacity);
    stars->Vz = safe_realloc(stars->Vz, sizeof(double) * stars->capacity);
    stars->mass = safe_realloc(stars->mass, sizeof(double) * stars->capacity);
}

int read_file(char *filename, Star *stars) {
    FILE *file = fopen(filename, "r");
    if (!file) return -1;

    char line[MAX_LINE];

    while (fgets(line, sizeof(line), file)) {
        // Ignorar líneas vacías o encabezados
        if (line[0] == 's' || line[0] == '\n') continue;

        // Eliminar el salto de línea final, si lo hay
        line[strcspn(line, "\r\n")] = '\0';

        char *tokens[9];
        int i = 0;

        // Separar por comas
        char *saveptr; //str_tok no es thread safe
        char *token = strtok_r(line, DELIMITER, &saveptr);
        while (token && i < 9) {
            tokens[i++] = token;
            token = strtok_r(NULL, DELIMITER, &saveptr);
        }

        // Si no hay suficientes columnas, descartamos
        if (i != 9) continue;

        // Validar campos importantes (parallax, mean_g, color)
        if (strcmp(tokens[3], "null") == 0 ||
            strcmp(tokens[7], "null") == 0 ||
            strcmp(tokens[8], "null") == 0) {
            continue;
        }

        // Asegurar capacidad antes de escribir
        if (stars->size >= stars->capacity) {
            stars->capacity += 10000;
            realloc_stars(stars);
        }

        int idx = stars->size++;

        stars->id[idx] = strtoul(tokens[0], NULL, 10);
        stars->ra[idx] = strtod(tokens[1], NULL);
        stars->dec[idx] = strtod(tokens[2], NULL);
        stars->parallax[idx] = strtod(tokens[3], NULL);
        stars->pmra[idx] = strtod(tokens[4], NULL);
        stars->pmdec[idx] = strtod(tokens[5], NULL);

        // radial_velocity puede ser null
        stars->radial_velocity[idx] = strcmp(tokens[6], "null") == 0 ? 0.0 : strtod(tokens[6], NULL);

        stars->mean_g[idx] = strtof(tokens[7], NULL);
        stars->color[idx] = strtof(tokens[8], NULL);
    }

    fclose(file);
    return 0;
}

void complete_data(Star *stars) {
    const double R[3][3] = {
        {-0.05487556, -0.87343709, -0.48383502},
        {+0.49410943, -0.44482963, +0.74698225},
        {-0.86766615, -0.19807637, +0.45598380}
    };

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

            double norm_C = sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
            double V_radial_sin = (C[0] * V[0] + C[1] * V[1] + C[2] * V[2]) / norm_C;

            double R_gal = sqrt(C[0] * C[0] + C[1] * C[1]);
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
}

int getstarsfromfile(char *dirname, Star *stars) {
    struct timeval start, end;
    struct dirent **filelist;



    // Obtener la lista de archivos
    int num_files = scandir(dirname, &filelist, NULL, alphasort);
    if (num_files < 0) {
        perror("No se pudo abrir el directorio");
        return -1;
    }

    gettimeofday(&start, NULL);

    // Lista de archivos válidos para procesar
    char **valid_files = malloc(num_files * sizeof(char *));
    int valid_count = 0;

    for (int i = 0; i < num_files; i++) {
        if (filelist[i]->d_name[0] == '.') {
            free(filelist[i]);
            continue;
        }
        char path[1000];
        sprintf(path, "%s/%s", dirname, filelist[i]->d_name);
        valid_files[valid_count++] = strdup(path);
        free(filelist[i]);
    }
    free(filelist);
    //Se preasigna una estimación de memoria para los datos en función del número de archivos
    stars->capacity = valid_count*500000;
    stars->size = 0;
    realloc_stars(stars);
    printf("Iniciando lectura de %d archivos\n", valid_count);
#pragma omp parallel
    {
        Star *temp = malloc(sizeof(Star));
        memset(temp, 0, sizeof(Star));
        temp->capacity = 700000;
        realloc_stars(temp);
        temp->size = 0;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < valid_count; i++) {
            read_file(valid_files[i], temp);
            complete_data(temp);
            int start_idx;

            // Reservamos espacio exacto solo si hay algo para copiar
            if (temp->size > 0) {
#pragma omp critical
                {
                    if (stars->size + temp->size > stars->capacity) {
                        stars->capacity = stars->size + temp->size + 1000000;
                        realloc_stars(stars);
                    }
                    start_idx = stars->size;
                    stars->size += temp->size;
                }

                // Copiamos fuera del critical
                for (int j = 0; j < temp->size; j++) {
                    int idx = start_idx + j;
                    stars->id[idx] = temp->id[j];
                    stars->ra[idx] = temp->ra[j];
                    stars->dec[idx] = temp->dec[j];
                    stars->parallax[idx] = temp->parallax[j];
                    stars->pmra[idx] = temp->pmra[j];
                    stars->pmdec[idx] = temp->pmdec[j];
                    stars->radial_velocity[idx] = temp->radial_velocity[j];
                    stars->mean_g[idx] = temp->mean_g[j];
                    stars->color[idx] = temp->color[j];
                    stars->Cx[idx] = temp->Cx[j];
                    stars->Cy[idx] = temp->Cy[j];
                    stars->Cz[idx] = temp->Cz[j];
                    stars->Vx[idx] = temp->Vx[j];
                    stars->Vy[idx] = temp->Vy[j];
                    stars->Vz[idx] = temp->Vz[j];
                    stars->mass[idx] = temp->mass[j];
                }
            }
            temp->size = 0; // limpio para siguiente archivo
        }
        free_stars(temp);
    }
    free(valid_files);
    gettimeofday(&end, NULL);

    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    printf("\nLeídas y trasladadas %d estrellas a memoria ocupando %.2f MB en %.2f segundos\n",
           stars->size,
           (stars->capacity * sizeof(double) * 11 + stars->capacity * sizeof(float) * 2 + stars->capacity * sizeof(
                unsigned long)) / (1024.0 * 1024.0), seconds);
    return stars->size;
}
