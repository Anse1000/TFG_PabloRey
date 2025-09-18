#ifndef AUX_FUN_H
#define AUX_FUN_H
#include "types.h"
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif
void *safe_realloc(void *ptr, size_t size);

void free_stars(Star *stars);

void resize_stars(Star *stars);

void free_aux(Star *stars);

double get_seconds(struct timeval start, struct timeval end);

static inline int get_octant(double cx, double cy, double cz, double x, double y, double z) {
    return ((x >= cx) << 2) | ((y >= cy) << 1) | (z >= cz);
}

void reorder_stars(Star *stars, double cx, double cy, double cz, unsigned int *offsets);

void compute_root_bounds(Star *stars, float *cx, float *cy, float *cz, float *hs, float *min_node_size,double min_subdivisions);

void write_results(Star *estrellas, const char *outputfile,const char*name,int step);

void estimate_dt(Star *stars, double *dt);
#ifdef __cplusplus
}
#endif
#endif //AUX_FUN_H
