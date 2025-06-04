#ifndef AUX_FUN_H
#define AUX_FUN_H
#include "types.h"

void free_stars(Star *stars);
void resize_stars(Star *stars);

void free_tree(Octree *tree);
void resize_tree(Octree *tree);

void free_aux(Star *stars);
double get_seconds(struct timeval start, struct timeval end);

#endif //AUX_FUN_H
