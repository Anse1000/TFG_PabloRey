#include "file_handler.h"
#include <dirent.h>
#include <omp.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "calculations.h"
#include "aux_fun.h"

void process_line(const char *line, Star *stars) {
    char *tokens[12];
    int i = 0;

    char *saveptr;
    char *token = strtok_r((char *)line, DELIMITER, &saveptr);
    while (token && i < 12) {
        tokens[i++] = token;
        token = strtok_r(NULL, DELIMITER, &saveptr);
    }

    // Validar que la línea tenga 12 columnas correctamente y que ciertos valores no sean "null"
    if (i == 12 && strcmp(tokens[11], "null") != 0 &&
        strcmp(tokens[7], "null") != 0 ){
        // Expandir arreglo si es necesario
        if (stars->size >= stars->capacity) {
            stars->capacity += 10000;
            resize_stars(stars);
        }

        unsigned long idx = stars->size++;

        stars->id[idx] = strtoul(tokens[0], NULL, 10);
        stars->ra[idx] = strtod(tokens[1], NULL);
        stars->dec[idx] = strtod(tokens[2], NULL);
        stars->pmra[idx] = strtod(tokens[3], NULL);
        stars->pmdec[idx] = strtod(tokens[4], NULL);
        stars->radial_velocity[idx] = strcmp(tokens[5], "null") == 0 ? 0.0 : strtod(tokens[5], NULL);
        stars->color[idx] = strtof(tokens[7], NULL);

        // Validar rango de color
        if (stars->color[idx] < 0.3 || stars->color[idx] > 2) {
            stars->size--;
            return;
        }

        stars->mass[idx] = strcmp(tokens[8], "null") == 0 ? 0.0F : strtof(tokens[8], NULL);
        stars->distance[idx] = strtod(tokens[11], NULL);
        }
}

int read_file(const char *filename, Star *stars) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error abriendo el archivo");
        return -1;
    }

    // Asignamos un búfer con suficiente tamaño para manejar fragmentos
    char *buffer = malloc(READ_BLOCK_SIZE + 1); // +1 para la terminación nula
    if (!buffer) {
        perror("Error asignando memoria al búfer");
        fclose(file);
        return -1;
    }

    size_t leftover = 0; // Bytes restantes (línea incompleta)
    size_t bytes_read; // Bytes leídos en cada iteración
    char *line_start;
    char *newline;

    while ((bytes_read = fread(buffer + leftover, 1, READ_BLOCK_SIZE - leftover, file)) > 0) {
        bytes_read += leftover; // Considerar el sobrante de la iteración anterior
        buffer[bytes_read] = '\0'; // Asegurarnos de que el búfer esté finalizado en cada lectura

        line_start = buffer; // Inicio de la línea actual

        // Buscar las líneas completas dentro del bloque leído
        while ((newline = strchr(line_start, '\n')) != NULL) {
            *newline = '\0'; // Finalizar línea actual

            // Procesar la línea si no es encabezado u hoja vacía
            if (line_start[0] != 's' && line_start[0] != '\0') {
                process_line(line_start, stars);
            }

            // Mover al siguiente inicio de línea
            line_start = newline + 1;
        }

        // Manejo del sobrante (línea cortada) al comienzo del búfer
        leftover = strlen(line_start);
        if (leftover > 0) {
            if (leftover > READ_BLOCK_SIZE) {
                fprintf(stderr, "Error: línea demasiado grande para el búfer\n");
                free(buffer);
                fclose(file);
                return -1;
            }
            memmove(buffer, line_start, leftover);
        }
    }

    // Procesar última línea si no termina en '\n'
    if (leftover > 0) {
        buffer[leftover] = '\0';
        if (buffer[0] != 's' && buffer[0] != '\0') {
            process_line(buffer, stars);
        }
    }

    free(buffer);
    fclose(file);
    return 0;
}

unsigned long getstarsfromfile(char *dirname, Star *stars) {
    struct timeval start, end;
    struct dirent **filelist;

    // Obtener la lista de archivos
    int num_files = scandir(dirname, &filelist, NULL, alphasort);
    if (num_files < 0) {
        perror("No se pudo abrir el directorio");
        return -1;
    }

    gettimeofday(&start, NULL);

    // Lista de archivos válidos para procesar
    char **valid_files = malloc(num_files * sizeof(char *));
    int valid_count = 0;

    for (int i = 0; i < num_files; i++) {
        if (filelist[i]->d_name[0] == '.') {
            free(filelist[i]);
            continue;
        }
        char path[1000];
        sprintf(path, "%s/%s", dirname, filelist[i]->d_name);
        valid_files[valid_count++] = strdup(path);
        free(filelist[i]);
    }
    free(filelist);
    //Se preasigna una estimación de memoria para los datos en función del número de archivos
    stars->capacity = valid_count * 400000;
    stars->size = 0;
    resize_stars(stars);
    printf("\nIniciando lectura de %d archivos usando %d threads\n", valid_count, omp_get_max_threads());
    fflush(stdout);
#pragma omp parallel
    {
        Star *temp = malloc(sizeof(Star));
        memset(temp, 0, sizeof(Star));
        temp->capacity = 700000;
        resize_stars(temp);
        temp->size = 0;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < valid_count; i++) {
            read_file(valid_files[i], temp);
            complete_data(temp);
            unsigned long start_idx;
            // Reservamos espacio exacto solo si hay algo para copiar
            if (temp->size > 0) {
#pragma omp critical
                {
                    if (stars->size + temp->size > stars->capacity) {
                        stars->capacity = stars->size + temp->size + 1000000;
                        resize_stars(stars);
                    }
                    start_idx = stars->size;
                    stars->size += temp->size;
                }

                // Copiamos fuera del critical
                for (size_t j = 0; j < temp->size; j++) {
                    size_t idx = start_idx + j;
                    stars->id[idx] = temp->id[j];
                    stars->ra[idx] = temp->ra[j];
                    stars->dec[idx] = temp->dec[j];
                    stars->distance[idx] = temp->distance[j];
                    stars->pmra[idx] = temp->pmra[j];
                    stars->pmdec[idx] = temp->pmdec[j];
                    stars->radial_velocity[idx] = temp->radial_velocity[j];
                    stars->color[idx] = temp->color[j];
                    stars->Cx[idx] = temp->Cx[j];
                    stars->Cy[idx] = temp->Cy[j];
                    stars->Cz[idx] = temp->Cz[j];
                    stars->Vx[idx] = temp->Vx[j];
                    stars->Vy[idx] = temp->Vy[j];
                    stars->Vz[idx] = temp->Vz[j];
                    stars->mass[idx] = temp->mass[j];
                }
            }
            temp->size = 0; // limpio para siguiente archivo
        }
        free_stars(temp);
    }
    stars->capacity = stars->size;
    resize_stars(stars);
    free(valid_files);
    gettimeofday(&end, NULL);

    double seconds = get_seconds(start, end);
    printf("Leídas y trasladadas %lu estrellas a memoria ocupando %.2f MB en %.2f segundos\n",
           stars->size,
           (stars->capacity * sizeof(double) * 13 + stars->capacity * sizeof(float) * 1 + stars->capacity * sizeof(
                unsigned long)) / (1024.0 * 1024.0), seconds);
    fflush(stdout);
    return stars->size;
}

void write_hdf5_chunks(Star *estrellas,
                       const char *directory,
                       const char *base_filename,
                       unsigned int num_chunks,
                       const size_t *chunk_sizes,
                       const size_t *chunk_offsets)
{
    struct stat st = {0};
    if (stat(directory, &st) == -1) mkdir(directory, 0755);

    #pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < num_chunks; i++)
    {
        size_t start = chunk_offsets[i];
        size_t count = chunk_sizes[i];

        char *filename = malloc(strlen(directory) + strlen(base_filename) + 20);
        sprintf(filename, "%s/%s_%02u.h5", directory, base_filename, i);

        hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_sec2(plist); // óptimo para Lustre
        hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
        H5Pclose(plist);

        hsize_t dims[1] = { count };
        hid_t space = H5Screate_simple(1, dims, NULL);

        #define WRITE_DATASET(name, type, ptr) \
            do { \
                hid_t dset = H5Dcreate(file_id, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
                H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (ptr)+start); \
                H5Dclose(dset); \
            } while(0)

        WRITE_DATASET("ID",   H5T_NATIVE_UINT64, estrellas->id);
        WRITE_DATASET("X",    H5T_NATIVE_DOUBLE, estrellas->Cx);
        WRITE_DATASET("Y",    H5T_NATIVE_DOUBLE, estrellas->Cy);
        WRITE_DATASET("Z",    H5T_NATIVE_DOUBLE, estrellas->Cz);
        WRITE_DATASET("MASS", H5T_NATIVE_FLOAT, estrellas->mass);

        H5Sclose(space);
        H5Fclose(file_id);
        free(filename);
    }
}

// ----------------------------
// Write master HDF5
// ----------------------------
void write_master_hdf5(const char *directory,
                       const char *base_filename,
                       int step,
                       unsigned int num_chunks,
                       const size_t *chunk_sizes,
                       const size_t *id_min,
                       const size_t *id_max)
{
    char *filename = malloc(strlen(directory) + strlen(base_filename) + 20);
    sprintf(filename, "%s/step_%04d_master.h5", directory, step+1);

    hid_t f = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // No necesitamos free(filename)

    hid_t g = H5Gcreate(f, "/global", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    size_t total = 0;
    for (unsigned int i = 0; i < num_chunks; i++) total += chunk_sizes[i];

    H5LTset_attribute_ulong(f, "/global", "total_stars", &total, 1);
    H5LTset_attribute_uint(f, "/global", "num_chunks", &num_chunks, 1);
    H5Gclose(g);

    hid_t gc = H5Gcreate(f, "/chunks", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (unsigned int i = 0; i < num_chunks; i++)
    {
        char path[64];
        snprintf(path, sizeof(path), "/chunks/chunk_%02u", i);
        hid_t gk = H5Gcreate(f, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        char chunkfile[1024];
        snprintf(chunkfile, sizeof(chunkfile), "%s_%02u.h5", base_filename, i);

        H5LTset_attribute_string(f, path, "filename", chunkfile);
        H5LTset_attribute_ulong(f, path, "size", &chunk_sizes[i], 1);
        H5LTset_attribute_ulong(f, path, "id_min", &id_min[i], 1);
        H5LTset_attribute_ulong(f, path, "id_max", &id_max[i], 1);

        // No necesitamos free(chunkfile)
        H5Gclose(gk);
    }

    H5Gclose(gc);
    H5Fclose(f);
    free(filename);
}

// ----------------------------
// Write master XDMF
// ----------------------------
void write_master_xdmf(const char *directory,
                       const char *base_filename,
                       int step,
                       unsigned int num_chunks,
                       const size_t *chunk_sizes)
{
    char *filename = malloc(strlen(directory) + strlen(base_filename) + 20);
    sprintf(filename,
             "%s/step_%04d.xmf",
             directory, step+1);

    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Error abriendo XDMF maestro");
        return;
    }

    fprintf(f,
"<?xml version=\"1.0\" ?>\n"
"<Xdmf Version=\"3.0\">\n"
"  <Domain>\n"
"    <Grid Name=\"Stars\" GridType=\"Collection\" CollectionType=\"Spatial\">\n");

    for (unsigned int i = 0; i < num_chunks; i++)
    {
        fprintf(f,
"      <Grid Name=\"chunk_%02u\" GridType=\"Uniform\">\n"
"        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%zu\"/>\n"
"        <Geometry GeometryType=\"XYZ\">\n"
"          <DataItem Dimensions=\"%zu\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
"            %s_%02u.h5:/X\n"
"          </DataItem>\n"
"          <DataItem Dimensions=\"%zu\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
"            %s_%02u.h5:/Y\n"
"          </DataItem>\n"
"          <DataItem Dimensions=\"%zu\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
"            %s_%02u.h5:/Z\n"
"          </DataItem>\n"
"        </Geometry>\n"
"        <Attribute Name=\"MASS\" AttributeType=\"Scalar\" Center=\"Node\">\n"
"          <DataItem Dimensions=\"%zu\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n"
"            %s_%02u.h5:/MASS\n"
"          </DataItem>\n"
"        </Attribute>\n"
"        <Attribute Name=\"ID\" AttributeType=\"Scalar\" Center=\"Node\">\n"
"          <DataItem Dimensions=\"%zu\" NumberType=\"UInt\" Precision=\"8\" Format=\"HDF\">\n"
"            %s_%02u.h5:/ID\n"
"          </DataItem>\n"
"        </Attribute>\n"
"      </Grid>\n",
            i,
            chunk_sizes[i],
            chunk_sizes[i], base_filename, i,
            chunk_sizes[i], base_filename, i,
            chunk_sizes[i], base_filename, i,
            chunk_sizes[i], base_filename, i,
            chunk_sizes[i], base_filename, i
        );
    }

    fprintf(f,
"    </Grid>\n"
"  </Domain>\n"
"</Xdmf>\n");

    fclose(f);
    free(filename);
}


// ----------------------------
// Main write_results
// ----------------------------
void write_results_hdf5(Star *estrellas,
                        const char *output_dir,
                        const char *name,
                        int step)
{
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) mkdir(output_dir, 0755);

    char *stepdir = malloc(strlen(output_dir) + strlen(name) + 10);
    sprintf(stepdir, "%s/step_%04d", output_dir, step+1);
    mkdir(stepdir, 0755);

    const unsigned int num_chunks = 50;
    size_t *chunk_sizes   = malloc(num_chunks * sizeof(size_t));
    size_t *chunk_offsets = malloc(num_chunks * sizeof(size_t));
    size_t *id_min        = malloc(num_chunks * sizeof(size_t));
    size_t *id_max        = malloc(num_chunks * sizeof(size_t));

    size_t N = estrellas->size;
    size_t base = N / num_chunks;
    size_t rem  = N % num_chunks;
    size_t offset = 0;

    for (unsigned int i = 0; i < num_chunks; i++)
    {
        chunk_sizes[i]   = base + (i < rem ? 1 : 0);
        chunk_offsets[i] = offset;
        id_min[i] = estrellas->id[offset];
        id_max[i] = estrellas->id[offset + chunk_sizes[i] - 1];
        offset += chunk_sizes[i];
    }

    write_hdf5_chunks(estrellas, stepdir, name,
                      num_chunks, chunk_sizes, chunk_offsets);
    write_master_hdf5(stepdir, name, step,
                      num_chunks, chunk_sizes, id_min, id_max);
    write_master_xdmf(stepdir, name, step,
                      num_chunks, chunk_sizes);

    free(chunk_sizes);
    free(chunk_offsets);
    free(id_min);
    free(id_max);
    free(stepdir);
}
