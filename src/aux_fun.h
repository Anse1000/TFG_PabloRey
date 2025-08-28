#ifndef AUX_FUN_H
#define AUX_FUN_H
#include "types.h"
#include <sys/time.h>
#include "octree.h"

#ifdef __cplusplus
extern "C" {
#endif
void free_stars(Star *stars);

void resize_stars(Star *stars);

void resize_tree(Octree *tree);

void free_aux(Star *stars);

void free_tree(Octree *tree);

double get_seconds(struct timeval start, struct timeval end);

static inline int get_octant(double cx, double cy, double cz, double x, double y, double z) {
    return ((x >= cx) << 2) | ((y >= cy) << 1) | (z >= cz);
}

void reorder_stars(Star *stars, double cx, double cy, double cz, unsigned int *offsets);

void compute_root_bounds(Star *stars, float *cx, float *cy, float *cz, float *hs, float *min_node_size,double min_subdivisions);

void write_chunks(Star *estrellas, const char *base_filename, const char *directory,unsigned int num_chunks,int add_mass);

void write_results(Star *estrellas, const char *outputfile,const char*name,int steps,int step);

void estimate_dt(Star *stars, double *dt);
#ifdef __cplusplus
}
#endif
#endif //AUX_FUN_H
