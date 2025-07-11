#include "simulation.h"
#include <omp.h>
#include "aux_fun.h"

double MIN_NODE_SIZE=1e-10;

//funcion de prueba: calcula la aceleracion de una sola estrella con TODAS
void compute_aceleration_single(const Star *stars, double *ax, double *ay, double *az, const unsigned long index,
                                double *seconds) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    for (unsigned long i = 0; i < stars->size; i++) {
        if (i != index) {
            //calcular distancia a la estrella
            double dx = stars->Cx[i] - stars->Cx[index];
            double dy = stars->Cy[i] - stars->Cy[index];
            double dz = stars->Cz[i] - stars->Cz[index];
            double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;
            double dist = sqrt(dist_sq);
            //calcular fuerza aplicada a la estrella
            double force = -G * stars->mass[i] / (dist_sq * dist);
            *ax = fma(force, dx, *ax);
            *ay = fma(force, dy, *ay);
            *az = fma(force, dz, *az);
        }
    }
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

long octree_new_node(Octree *tree, float cx, float cy, float cz, float half_size) {
    if (tree->size >= tree->capacity) {
        tree->capacity *= 1.2;
        resize_tree(tree);
    }
    size_t i = tree->size++;

    tree->center_x[i] = cx;
    tree->center_y[i] = cy;
    tree->center_z[i] = cz;
    tree->half_size[i] = half_size;

    tree->mass[i] = 0.0F;
    tree->com_x[i] = 0.0;
    tree->com_y[i] = 0.0;
    tree->com_z[i] = 0.0;
    tree->star_index[i] = -1;

    for (int j = 0; j < 8; j++)
        tree->children[i][j] = INVALID_INDEX;
    return i;
}

inline int get_octant(double cx, double cy, double cz, double x, double y, double z) {
    return ((x >= cx) << 2) | ((y >= cy) << 1) | (z >= cz);
}

void octree_insert(Octree *tree, Star *stars, long node_index, long star_index) {
    float cx = tree->center_x[node_index];
    float cy = tree->center_y[node_index];
    float cz = tree->center_z[node_index];
    float hs = tree->half_size[node_index];

    double x = stars->Cx[star_index];
    double y = stars->Cy[star_index];
    double z = stars->Cz[star_index];
    float m = stars->mass[star_index];

    // Actualizar masa y centro de masa del nodo
    float old_mass = tree->mass[node_index];
    float new_mass = old_mass + m;

    tree->com_x[node_index] = (tree->com_x[node_index] * old_mass + x * m) / new_mass;
    tree->com_y[node_index] = (tree->com_y[node_index] * old_mass + y * m) / new_mass;
    tree->com_z[node_index] = (tree->com_z[node_index] * old_mass + z * m) / new_mass;
    tree->mass[node_index] = new_mass;

    // Si el nodo es demasiado pequeño, no subdividir más
    if (hs * 2.0 <= MIN_NODE_SIZE) {
        if (tree->star_index[node_index] == -1)
            tree->star_index[node_index] = star_index;  // asignar primera estrella
        // si ya hay una, se quedan varias aquí (no se subdivide más)
        return;
    }

    int oct = get_octant(cx, cy, cz, x, y, z);

    if (tree->children[node_index][oct] == INVALID_INDEX) {
        // Crear nuevo nodo hijo
        float offset = hs * 0.5F;
        float new_cx = cx + ((oct & 4) ? offset : -offset);
        float new_cy = cy + ((oct & 2) ? offset : -offset);
        float new_cz = cz + ((oct & 1) ? offset : -offset);

        long child_index = octree_new_node(tree, new_cx, new_cy, new_cz, offset);
        tree->children[node_index][oct] = child_index;

        // Insertar directamente en el hijo
        tree->star_index[child_index] = star_index;
        tree->mass[child_index] = m;
        tree->com_x[child_index] = x;
        tree->com_y[child_index] = y;
        tree->com_z[child_index] = z;
    } else {
        long child = tree->children[node_index][oct];
        if (tree->star_index[child] >= 0) {
            long existing_star = tree->star_index[child];
            tree->star_index[child] = -1;

            // Resetear propiedades acumuladas del hijo antes de reinserciones
            tree->mass[child] = 0.0F;
            tree->com_x[child] = 0.0;
            tree->com_y[child] = 0.0;
            tree->com_z[child] = 0.0;

            octree_insert(tree, stars, child, existing_star);
            octree_insert(tree, stars, child, star_index);
        } else {
            octree_insert(tree, stars, child, star_index);
        }
    }
}

void compute_acceleration_bh(const Star *stars, const Octree *tree,
                             long node_idx, long star_idx, double theta,
                             double *ax, double *ay, double *az) {
    // Ignorar si es una hoja con la misma estrella
    if (tree->star_index[node_idx] == star_idx)
        return;

    // Diferencia de posición
    double dx = tree->com_x[node_idx] - stars->Cx[star_idx];
    double dy = tree->com_y[node_idx] - stars->Cy[star_idx];
    double dz = tree->com_z[node_idx] - stars->Cz[star_idx];
    double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;
    double dist = sqrt(dist_sq);

    double s = 2.0 * tree->half_size[node_idx]; // ancho total del nodo

    if ((s / dist) < theta || tree->star_index[node_idx] >= 0) {
        // Tratar como nodo lejano o hoja
        double force = -G * tree->mass[node_idx] / (dist_sq * dist);
        *ax = fma(force, dx, *ax);
        *ay = fma(force, dy, *ay);
        *az = fma(force, dz, *az);
    } else {
        // Recursión en hijos
        for (int i = 0; i < 8; i++) {
            long child = tree->children[node_idx][i];
            if (child != INVALID_INDEX) {
                compute_acceleration_bh(stars, tree, child, star_idx, theta, ax, ay, az);
            }
        }
    }
}

void aux_time_bh(const Star *stars, const Octree *tree, long node_idx, long index, double theta, double *ax,
                 double *ay, double *az, double *seconds)  {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    compute_acceleration_bh(stars, tree, node_idx, index, theta, ax, ay, az);
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

Octree *build_tree(Star *stars) {
    struct timeval start, end;
    size_t initial_capacity = 10000;
    gettimeofday(&start, NULL);
    printf("Iniciando construccion del Arbol\n");
    fflush(stdout);

    Octree *tree = malloc(sizeof(Octree));
    memset(tree, 0, sizeof(Octree));

    tree->capacity = initial_capacity;
    tree->size = 0;
    resize_tree(tree);

    for (size_t i = 0; i < initial_capacity; i++) {
        for (int j = 0; j < 8; j++) tree->children[i][j] = INVALID_INDEX;
        tree->star_index[i] = -1;
    }

    float cx, cy, cz;
    float hs;
    compute_root_bounds(stars, &cx, &cy, &cz, &hs,&MIN_NODE_SIZE);
    printf("MIN_NODE_SIZE = %f\n", MIN_NODE_SIZE); fflush(stdout);

    long root = octree_new_node(tree, cx, cy, cz, hs);

    for (unsigned long i = 0; i < stars->size; i++) {
        octree_insert(tree, stars, root, i);
    }
    if (tree->size < tree->capacity) {
        tree->capacity = tree->size;
        resize_tree(tree);
    }
    gettimeofday(&end, NULL);
    size_t memory = tree->capacity * (
                        sizeof(double) * 3 + // center_x, center_y, center_z, half_size, com_x, com_y, com_z
                        sizeof(double) + // mass
                        sizeof(float) * 4 +
                        sizeof(unsigned int[8]) + // children (8 longs por nodo)
                        sizeof(long) // star_index
                    );
    double secs = get_seconds(start, end);
    printf("Árbol de %ld nodos creado en %.4f segundos usando %lu MB \n", tree->capacity, secs, memory / 1024 / 1024);
    fflush(stdout);
    return tree;
}

// Función principal de simulación
void simulate(Star *estrellas, const long N, const char* outputfile) {
    struct timeval start, end;
    double DT2 = 0.5 * DT;
    double *ax = malloc(N * sizeof(double));
    double *ay = malloc(N * sizeof(double));
    double *az = malloc(N * sizeof(double));

    gettimeofday(&start, NULL);
    for (int step = 0; step < STEPS; step++) {
        Octree *octree = build_tree(estrellas);
        printf("Iniciando fase 1\n"); fflush(stdout);
        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            compute_acceleration_bh(estrellas, octree, 0, i, 0.2, &ax[i], &ay[i], &az[i]);
            // Leapfrog integration: actualizar velocidad a mitad de paso
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
            // Actualizar posición
            estrellas->Cx[i] = fma(DT, estrellas->Vx[i], estrellas->Cx[i]);
            estrellas->Cy[i] = fma(DT, estrellas->Vy[i], estrellas->Cy[i]);
            estrellas->Cz[i] = fma(DT, estrellas->Vz[i], estrellas->Cz[i]);
        }
        free_tree(octree);
        //Reconstruir con nuevas posiciones
        octree = build_tree(estrellas);
        printf("Iniciando fase 2\n"); fflush(stdout);
        #pragma omp parallel for
        for (long i = 0; i < N; i++) {
            compute_acceleration_bh(estrellas, octree, 0, i, 0.2, &ax[i], &ay[i], &az[i]);
            // Completar actualización de velocidad
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
        }
        printf("Paso %d realizado\n", step + 1);
        fflush(stdout);
    }
    gettimeofday(&end, NULL);
    double seconds = get_seconds(start, end);
    int hours = (int)(seconds / 3600);
    int minutes = ((int)seconds % 3600) / 60;
    double remaining_seconds = fmod(seconds, 60.0);
    printf("Simuladas %ld estrellas en %02d:%02d:%05.2f (hh:mm:ss) usando %d threads\n",
       N, hours, minutes, remaining_seconds, omp_get_max_threads());
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
    FILE *file = fopen(outputfile, "w");
    if (!file) {
        fprintf(stderr, "Error al abrir el archivo de salida\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n", estrellas->id[i], estrellas->Cx[i], estrellas->Cy[i],
                estrellas->Cz[i]);
    }
}