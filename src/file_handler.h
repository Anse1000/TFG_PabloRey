#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H
#include "types.h"
#include <dirent.h>

#define MAX_LINE 2000
#define DELIMITER ","


int getstarsfromfile(char *dirname, Star *estrellas);

#endif // FILE_HANDLER_H
