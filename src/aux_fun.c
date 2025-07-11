#include "aux_fun.h"

static void *safe_realloc(void *ptr, const size_t size) {
    void *tmp = realloc(ptr, size);
    if (!tmp) {
        perror("Fallo al hacer el resize de memoria");
        exit(EXIT_FAILURE);
    }
    return tmp;
}
double get_seconds(const struct timeval start, const struct timeval end) {
    time_t diff = end.tv_sec - start.tv_sec;
    suseconds_t diff_us = end.tv_usec - start.tv_usec;
    return diff + diff_us / 1000000.0;
}

void free_stars(Star *stars) {
    free(stars->id);
    free(stars->Cx);
    free(stars->Cy);
    free(stars->Cz);
    free(stars->Vx);
    free(stars->Vy);
    free(stars->Vz);
    free(stars->mass);
    free(stars);
}

void resize_stars(Star *stars) {
    stars->id = safe_realloc(stars->id, sizeof(long) * stars->capacity);
    stars->ra = safe_realloc(stars->ra, sizeof(double) * stars->capacity);
    stars->dec = safe_realloc(stars->dec, sizeof(double) * stars->capacity);
    stars->distance = safe_realloc(stars->distance, sizeof(double) * stars->capacity);
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
    stars->radius = safe_realloc(stars->radius, sizeof(float) * stars->capacity);
    stars->gravity = safe_realloc(stars->gravity, sizeof(float) * stars->capacity);
}
void resize_tree(Octree *tree) {
    tree->center_x = safe_realloc(tree->center_x, sizeof(float) * tree->capacity);
    tree->center_y = safe_realloc(tree->center_y, sizeof(float) * tree->capacity);
    tree->center_z = safe_realloc(tree->center_z, sizeof(float) * tree->capacity);
    tree->half_size = safe_realloc(tree->half_size, sizeof(float) * tree->capacity);
    tree->mass = safe_realloc(tree->mass, sizeof(float) * tree->capacity);
    tree->com_x = safe_realloc(tree->com_x, sizeof(double) * tree->capacity);
    tree->com_y = safe_realloc(tree->com_y, sizeof(double) * tree->capacity);
    tree->com_z = safe_realloc(tree->com_z, sizeof(double) * tree->capacity);
    tree->children = safe_realloc(tree->children, sizeof(long[8]) * tree->capacity);
    tree->star_index = safe_realloc(tree->star_index, sizeof(long) * tree->capacity);
}
void compute_root_bounds(Star *estrellas, float *center_x, float *center_y, float *center_z, float *half_size,double *min_node_size) {
    // Inicializar límites
    double min_cx = estrellas->Cx[0], max_cx = estrellas->Cx[0];
    double min_cy = estrellas->Cy[0], max_cy = estrellas->Cy[0];
    double min_cz = estrellas->Cz[0], max_cz = estrellas->Cz[0];

    // Calcular límites de posición
    for (unsigned long i = 0; i < estrellas->size; i++) {
        if (estrellas->Cx[i] < min_cx) min_cx = estrellas->Cx[i];
        else if (estrellas->Cx[i] > max_cx) max_cx = estrellas->Cx[i];
        if (estrellas->Cy[i] < min_cy) min_cy = estrellas->Cy[i];
        else if (estrellas->Cy[i] > max_cy) max_cy = estrellas->Cy[i];
        if (estrellas->Cz[i] < min_cz) min_cz = estrellas->Cz[i];
        else if (estrellas->Cz[i] > max_cz) max_cz = estrellas->Cz[i];
    }

    // Calcular centro
    *center_x = 0.5F * (min_cx + max_cx);
    *center_y = 0.5F * (min_cy + max_cy);
    *center_z = 0.5F * (min_cz + max_cz);

    // Calcular rango máximo
    double dx = max_cx - min_cx;
    double dy = max_cy - min_cy;
    double dz = max_cz - min_cz;
    double max_range = fmax(dx, fmax(dy, dz));

    // Usar margen de seguridad (20%) y dividir entre 2
    *half_size = 0.5F * max_range * 1.2F;

    //Elegir precisión para subdivisiones
    *min_node_size = max_range * 1e-7;
}
void free_tree(Octree *tree) {
    if (!tree) return;
    free(tree->center_x);
    free(tree->center_y);
    free(tree->center_z);
    free(tree->half_size);
    free(tree->mass);
    free(tree->com_x);
    free(tree->com_y);
    free(tree->com_z);
    free(tree->children);
    free(tree->star_index);
}

void free_aux(Star *estrellas) {
    free(estrellas->ra);
    free(estrellas->dec);
    free(estrellas->pmdec);
    free(estrellas->pmra);
    free(estrellas->color);
    free(estrellas->radius);
    free(estrellas->radial_velocity);
    free(estrellas->distance);
    free(estrellas->gravity);
    free(estrellas->mean_g);
    size_t total_bytes = estrellas->size * (
                             sizeof(double) * 6 + // ra, dec, pmdec, pmra, radial_velocity, distance
                             sizeof(float) * 4 // color, radius, gravity, mean_g
                         );
    printf("Liberados %.2lu MB de recursos auxiliares\n",total_bytes/1024/1024);
}
void swap_star_elements(Star *star, size_t i, size_t j) {
#define SWAP(arr) do { typeof((arr)[0]) tmp = (arr)[i]; (arr)[i] = (arr)[j]; (arr)[j] = tmp; } while (0)

    SWAP(star->id);

    SWAP(star->Cx);
    SWAP(star->Cy);
    SWAP(star->Cz);

    SWAP(star->Vx);
    SWAP(star->Vy);
    SWAP(star->Vz);

    SWAP(star->mass);

#undef SWAP
}