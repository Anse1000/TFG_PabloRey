#ifndef AUX_FUN_H
#define AUX_FUN_H
#include "types.h"

void free_stars(Star *stars);
void resize_stars(Star *stars);

void free_tree(Octree *tree);
void resize_tree(Octree *tree);

void free_aux(Star *stars);
#ifdef __cplusplus
extern "C" {
#endif
double get_seconds(struct timeval start, struct timeval end);

void compute_root_bounds(Star *estrellas, float *center_x, float *center_y, float *center_z, float *half_size,double *min_node_size);
void swap_star_elements(Star *star, size_t i, size_t j);
#ifdef __cplusplus
}
#endif
#endif //AUX_FUN_H
