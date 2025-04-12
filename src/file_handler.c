#include "file_handler.h"

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
            *allocated_size += 10000; // Aumentar de a 10000 elementos
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
            {&current->id, TYPE_UL},
            {&current->ra, TYPE_D},
            {&current->dec, TYPE_D},
            {&current->parallax, TYPE_D},
            {&current->pmra, TYPE_D},
            {&current->pmdec, TYPE_D},
            {&current->radial_velocity, TYPE_D},
            {&current->mean_g, TYPE_F},
            {&current->color, TYPE_F}
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
int getstarsfromfile(char *dirname, Star **estrellas) {
    struct timeval start, end;
    char path[1000];
    struct dirent *filedir;
    int dirname_length = strlen(dirname);
    int N = 0, count = 0, allocated = 0;
    int status = 0;
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
        status = read_file(path, estrellas, &N, &allocated);
        if (status != 0) {
            printf("Error en la lectura de archivos\n");
            return -1;
        }
        count++;
        if (count % 10 == 0) {
            printf("Leidos %d archivos...\n", count);
            fflush(stdout);
        }
    }
    //Ajustar memoria a tamaño exacto
    Star *temp = realloc(*estrellas, N * sizeof(Star)); // Redimensionar al número exacto de estrellas leídas
    if (temp != NULL) {
        *estrellas = temp;
    } else {
        perror("No se pudo redimensionar la memoria.\n");
    }
    gettimeofday(&end, NULL);
    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("\nLeídas y trasladadas %d estrellas a memoria ocupando %.2f MB en %.2f segundos\n", N,
           N * sizeof(Star) / (1024.0 * 1024.0), seconds);
    fflush(stdout);
    fflush(0);
    return N;
}