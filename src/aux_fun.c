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
    free(stars->ra);
    free(stars->dec);
    free(stars->distance);
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
    free(stars->radius);
    free(stars->gravity);
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
    tree->cx = safe_realloc(tree->cx, sizeof(double) * tree->capacity);
    tree->cy = safe_realloc(tree->cy, sizeof(double) * tree->capacity);
    tree->cz = safe_realloc(tree->cz, sizeof(double) * tree->capacity);
    tree->half_size = safe_realloc(tree->half_size, sizeof(float) * tree->capacity);
    tree->mass = safe_realloc(tree->mass, sizeof(float) * tree->capacity);
    tree->com_x = safe_realloc(tree->com_x, sizeof(float) * tree->capacity);
    tree->com_y = safe_realloc(tree->com_y, sizeof(float) * tree->capacity);
    tree->com_z = safe_realloc(tree->com_z, sizeof(float) * tree->capacity);
    tree->star_index = safe_realloc(tree->star_index, sizeof(long) * tree->capacity);
    tree->first_child_index = safe_realloc(tree->first_child_index, sizeof(long) * tree->capacity);

}
void free_tree(Octree *tree) {
    if (!tree) return;
    free(tree->cx);
    free(tree->cy);
    free(tree->cz);
    free(tree->half_size);
    free(tree->mass);
    free(tree->com_x);
    free(tree->com_y);
    free(tree->com_z);
    free(tree->star_index);
    free(tree->first_child_index);
    free(tree);
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
    size_t total_bytes = estrellas->size * (
                             sizeof(double) * 6 + // ra, dec, pmdec, pmra, radial_velocity, distance
                             sizeof(float) * 3 // color, radius, gravity
                         );
    printf("Liberados %.2lu MB de recursos auxiliares\n",total_bytes/1024/1024);
}