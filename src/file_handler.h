#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H
#include "types.h"
#include <dirent.h>
#include <omp.h>

#define MAX_LINE 2000
#define DELIMITER ","

#define TYPE_UL  0  // unsigned long
#define TYPE_D   1  // double
#define TYPE_F   2  // float

typedef struct {
    void *ptr;
    int type;
} FieldMap;

int getstarsfromfile(char *dirname, Star *estrellas);
void free_stars(Star *stars);

#endif // FILE_HANDLER_H
