#include "file_handler.h"

static void *safe_realloc(void *ptr, size_t size) {
    void *tmp = realloc(ptr, size);
    if (!tmp) {
        perror("Fallo al hacer el resize de memoria");
        exit(EXIT_FAILURE);
    }
    return tmp;
}

void realloc_stars(Star *stars) {
    stars->id = safe_realloc(stars->id,sizeof(long) * stars->capacity);
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
        if (line[0] == 's') continue;

        char *tokens[21];
        char *token = strtok(line, DELIMITER);
        int is_valid = 1;
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

        if (stars->size >= stars->capacity) {
            stars->capacity += 10000;
            realloc_stars(stars);
        }

        int idx = stars->size++;

        stars->id[idx]               = strtoul(tokens[0], NULL, 10);
        stars->ra[idx]               = strtod(tokens[1], NULL);
        stars->dec[idx]              = strtod(tokens[2], NULL);
        stars->parallax[idx]         = strtod(tokens[3], NULL);
        stars->pmra[idx]             = strtod(tokens[4], NULL);
        stars->pmdec[idx]            = strtod(tokens[5], NULL);
        stars->radial_velocity[idx]  = strtod(tokens[6], NULL);
        stars->mean_g[idx]           = strtof(tokens[7], NULL);
        stars->color[idx]            = strtof(tokens[8], NULL);
    }

    fclose(file);
    return 0;
}

int getstarsfromfile(char *dirname, Star *stars) {
    struct timeval start, end;
    char path[1000];
    struct dirent *filedir;
    int count = 0, status = 0;

    stars->size = 0;
    stars->capacity = 0;

    DIR *dir = opendir(dirname);
    if (dir == NULL) {
        perror("No se pudo abrir el directorio");
        return -1;
    }

    int dirname_length = strlen(dirname);
    strcpy(path, dirname);
    strcat(path, "/");

    gettimeofday(&start, NULL);

    while ((filedir = readdir(dir)) != NULL) {
        if (filedir->d_name[0] == '.') continue;

        strcpy(path + dirname_length + 1, filedir->d_name);

        status = read_file(path, stars);
        if (status != 0) {
            printf("Error en la lectura de archivos\n");
            closedir(dir);
            return -1;
        }

        count++;
        if (count % 10 == 0) {
            printf("Leídos %d archivos...\n", count);
            fflush(stdout);
        }
    }

    closedir(dir);
    gettimeofday(&end, NULL);

    double seconds = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    printf("\nLeídas y trasladadas %d estrellas a memoria ocupando %.2f MB en %.2f segundos\n",
           stars->size,
           (stars->capacity * sizeof(double) * 11 + stars->capacity * sizeof(float) * 2 + stars->capacity * sizeof(unsigned long)) / (1024.0 * 1024.0),
           seconds);
    return stars->size;
}
