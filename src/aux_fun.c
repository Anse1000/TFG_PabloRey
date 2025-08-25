#include "aux_fun.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "octree.h"

// Función mejorada para nodo raíz Barnes-Hut
void compute_root_bounds(Star *stars, float *cx, float *cy, float *cz, float *hs, float *min_node_size,double min_subdivisions) {
    if (stars->size == 0) return;

    // 1. Centro de masa inicial
    double total_mass = 0.0;
    double center_x = 0.0, center_y = 0.0, center_z = 0.0;

    for (size_t i = 0; i < stars->size; i++) {
        double m = stars->mass[i];
        total_mass += m;
        center_x += stars->Cx[i] * m;
        center_y += stars->Cy[i] * m;
        center_z += stars->Cz[i] * m;
    }

    center_x /= total_mass;
    center_y /= total_mass;
    center_z /= total_mass;

    // 2. Calcular bounding sphere inicial desde COM
    double max_dist_sq = 0.0;
    for (size_t i = 0; i < stars->size; i++) {
        double dx = stars->Cx[i] - center_x;
        double dy = stars->Cy[i] - center_y;
        double dz = stars->Cz[i] - center_z;
        double dist_sq = dx*dx + dy*dy + dz*dz;
        if (dist_sq > max_dist_sq) {
            max_dist_sq = dist_sq;
        }
    }

    double radius = sqrt(max_dist_sq);

    // 3. Algoritmo tipo Ritter: expandir esfera hacia puntos externos
    for (size_t i = 0; i < stars->size; i++) {
        double dx = stars->Cx[i] - center_x;
        double dy = stars->Cy[i] - center_y;
        double dz = stars->Cz[i] - center_z;
        double dist = sqrt(dx*dx + dy*dy + dz*dz);

        if (dist > radius) {
            // Mover el centro hacia este punto para mantener esfera compacta
            double shift = (dist - radius) * 0.5 / dist;
            center_x += dx * shift;
            center_y += dy * shift;
            center_z += dz * shift;
            radius += (dist - radius) * 0.5;
        }
    }

    // 4. Ajustar valores de salida
    *cx = (float)center_x;
    *cy = (float)center_y;
    *cz = (float)center_z;
    *hs = (float)(radius * 1.1); // margen 10%
    *min_node_size = *hs / min_subdivisions;
}

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
    tree->children = safe_realloc(tree->children, sizeof(unsigned int[8]) * tree->capacity);
    tree->star_index = safe_realloc(tree->star_index, sizeof(long) * tree->capacity);
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
//reordenar estrellas por octante para hacer calculos en gpu
void reorder_stars(Star *stars,double cx, double cy, double cz, unsigned int *offsets) {
    size_t counts[8] = {0};
    for (size_t i = 0; i < stars->size; i++) {
        int oct = get_octant(cx, cy, cz, stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        counts[oct]++;
    }
    offsets[0] = 0;
    for (int i = 1; i < 8; i++) {
        offsets[i] = offsets[i - 1] + counts[i - 1];
    }
    unsigned int ends[8];
    memcpy(ends, offsets, 8*sizeof(unsigned int));
    for (size_t i = 0; i < stars->size;) {
        int oct = get_octant(cx, cy, cz, stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        if (i >= offsets[oct] && i < ends[oct]) {
            // Ya está en su rango
            i++;
        } else {
            // Debe ir en ends[oct]
            size_t dest = ends[oct];
            swap_star_elements(stars, i, dest);
            ends[oct]++;
        }
    }
    printf("Estrellas ordenadas por octante\n"); fflush(stdout);
}