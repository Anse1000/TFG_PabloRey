#include "file_handler.h"
#include <dirent.h>
#include <omp.h>
#include "calculations.h"
#include "aux_fun.h"

void process_line(const char *line, Star *stars) {
    char *tokens[12];
    int i = 0;

    char *saveptr;
    char *token = strtok_r((char *)line, DELIMITER, &saveptr);
    while (token && i < 12) {
        tokens[i++] = token;
        token = strtok_r(NULL, DELIMITER, &saveptr);
    }

    // Validar que la línea tenga 12 columnas correctamente y que ciertos valores no sean "null"
    if (i == 12 && strcmp(tokens[11], "null") != 0 &&
        strcmp(tokens[7], "null") != 0 &&
        strcmp(tokens[6], "null") != 0) {
        // Expandir arreglo si es necesario
        if (stars->size >= stars->capacity) {
            stars->capacity += 10000;
            resize_stars(stars);
        }

        unsigned long idx = stars->size++;

        stars->id[idx] = strtoul(tokens[0], NULL, 10);
        stars->ra[idx] = strtod(tokens[1], NULL);
        stars->dec[idx] = strtod(tokens[2], NULL);
        stars->pmra[idx] = strtod(tokens[3], NULL);
        stars->pmdec[idx] = strtod(tokens[4], NULL);
        stars->radial_velocity[idx] = strcmp(tokens[5], "null") == 0 ? 0.0 : strtod(tokens[5], NULL);
        stars->mean_g[idx] = strtof(tokens[6], NULL);
        stars->color[idx] = strtof(tokens[7], NULL);

        // Validar rango de color
        if (stars->color[idx] < 0.3 || stars->color[idx] > 2) {
            stars->size--;
            return;
        }

        stars->mass[idx] = strcmp(tokens[8], "null") == 0 ? 0.0F : strtof(tokens[8], NULL);
        stars->radius[idx] = strcmp(tokens[9], "null") == 0 ? 0.0F : strtof(tokens[9], NULL);
        stars->gravity[idx] = strcmp(tokens[10], "null") == 0 ? 0.0F : strtof(tokens[10], NULL);
        stars->distance[idx] = strtod(tokens[11], NULL);
        }
}

int read_file(const char *filename, Star *stars) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error abriendo el archivo");
        return -1;
    }

    // Asignamos un búfer con suficiente tamaño para manejar fragmentos
    char *buffer = malloc(BLOCK_SIZE + 1); // +1 para la terminación nula
    if (!buffer) {
        perror("Error asignando memoria al búfer");
        fclose(file);
        return -1;
    }

    size_t leftover = 0; // Bytes restantes (línea incompleta)
    size_t bytes_read; // Bytes leídos en cada iteración
    char *line_start;
    char *newline;

    while ((bytes_read = fread(buffer + leftover, 1, BLOCK_SIZE - leftover, file)) > 0) {
        bytes_read += leftover; // Considerar el sobrante de la iteración anterior
        buffer[bytes_read] = '\0'; // Asegurarnos de que el búfer esté finalizado en cada lectura

        line_start = buffer; // Inicio de la línea actual

        // Buscar las líneas completas dentro del bloque leído
        while ((newline = strchr(line_start, '\n')) != NULL) {
            *newline = '\0'; // Finalizar línea actual

            // Procesar la línea si no es encabezado u hoja vacía
            if (line_start[0] != 's' && line_start[0] != '\0') {
                process_line(line_start, stars);
            }

            // Mover al siguiente inicio de línea
            line_start = newline + 1;
        }

        // Manejo del sobrante (línea cortada) al comienzo del búfer
        leftover = strlen(line_start);
        if (leftover > 0) {
            if (leftover > BLOCK_SIZE) {
                fprintf(stderr, "Error: línea demasiado grande para el búfer\n");
                free(buffer);
                fclose(file);
                return -1;
            }
            memmove(buffer, line_start, leftover);
        }
    }

    // Procesar última línea si no termina en '\n'
    if (leftover > 0) {
        buffer[leftover] = '\0';
        if (buffer[0] != 's' && buffer[0] != '\0') {
            process_line(buffer, stars);
        }
    }

    free(buffer);
    fclose(file);
    return 0;
}

unsigned long getstarsfromfile(char *dirname, Star *stars) {
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
    stars->capacity = valid_count * 400000;
    stars->size = 0;
    resize_stars(stars);
    printf("\nIniciando lectura de %d archivos usando %d threads\n", valid_count, omp_get_max_threads());
    fflush(stdout);
#pragma omp parallel
    {
        Star *temp = malloc(sizeof(Star));
        memset(temp, 0, sizeof(Star));
        temp->capacity = 700000;
        resize_stars(temp);
        temp->size = 0;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < valid_count; i++) {
            read_file(valid_files[i], temp);
            complete_data(temp);
            unsigned long start_idx;
            // Reservamos espacio exacto solo si hay algo para copiar
            if (temp->size > 0) {
#pragma omp critical
                {
                    if (stars->size + temp->size > stars->capacity) {
                        stars->capacity = stars->size + temp->size + 1000000;
                        resize_stars(stars);
                    }
                    start_idx = stars->size;
                    stars->size += temp->size;
                }

                // Copiamos fuera del critical
                for (size_t j = 0; j < temp->size; j++) {
                    size_t idx = start_idx + j;
                    stars->id[idx] = temp->id[j];
                    stars->ra[idx] = temp->ra[j];
                    stars->dec[idx] = temp->dec[j];
                    stars->distance[idx] = temp->distance[j];
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
                    stars->gravity[idx] = temp->gravity[j];
                    stars->radius[idx] = temp->radius[j];
                }
            }
            temp->size = 0; // limpio para siguiente archivo
        }
        free_stars(temp);
    }
    stars->capacity = stars->size;
    resize_stars(stars);
    free(valid_files);
    gettimeofday(&end, NULL);

    double seconds = get_seconds(start, end);
    printf("Leídas y trasladadas %lu estrellas a memoria ocupando %.2f MB en %.2f segundos\n",
           stars->size,
           (stars->capacity * sizeof(double) * 13 + stars->capacity * sizeof(float) * 4 + stars->capacity * sizeof(
                unsigned long)) / (1024.0 * 1024.0), seconds);
    fflush(stdout);
    return stars->size;
}
