#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H
#include "types.h"

#define DELIMITER ","
#define BLOCK_SIZE (4 * 1024 * 1024)  // 4MB

unsigned long getstarsfromfile(char *dirname, Star *estrellas);
void free_stars(Star *stars);

#endif // FILE_HANDLER_H
