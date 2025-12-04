#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H
#include "types.h"

#define DELIMITER ","
#define READ_BLOCK_SIZE (4 * 1024 * 1024)  // 4MB

unsigned long getstarsfromfile(char *dirname, Star *estrellas);
void free_stars(Star *stars);
#ifdef __cplusplus
extern "C" {
#endif
void write_results_hdf5(Star *estrellas, const char *output_dir, const char *name, int step);
#ifdef __cplusplus
}
#endif
#endif // FILE_HANDLER_H
