#include "simulation.h"
#include "aux_fun.h"

double MIN_NODE_SIZE=1e-10;

//Funcion de cómputo de aceleración
void compute_aceleration(Star *estrellas, double *ax, double *ay, double *az, int N) {
    for (int i = 0; i < N; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                //calcular distancia a la estrella
                double dx = estrellas->Cx[i] - estrellas->Cx[j];
                double dy = estrellas->Cy[i] - estrellas->Cy[j];
                double dz = estrellas->Cz[i] - estrellas->Cz[j];
                double dist_sq = dx * dx + dy * dy + dz * dz;
                double dist = sqrt(dist_sq);
                //calcular fuerza aplicada a la estrella
                double force = -G * estrellas->mass[j] / (dist_sq * dist);
                ax[i] = fma(force, dx, ax[i]);
                ay[i] = fma(force, dy, ay[i]);
                az[i] = fma(force, dz, az[i]);
            }
        }
    }
}

#ifdef AVX_512
void compute_aceleration_avx512(const Star *stars, double *ax, double *ay, double *az, int N) {
    const __m512d G_vec = _mm512_set1_pd(-G);

    for (int i = 0; i < N; i++) {
        __m512d xi = _mm512_set1_pd(stars->Cx[i]);
        __m512d yi = _mm512_set1_pd(stars->Cy[i]);
        __m512d zi = _mm512_set1_pd(stars->Cz[i]);

        __m512d sum_ax = _mm512_setzero_pd();
        __m512d sum_ay = _mm512_setzero_pd();
        __m512d sum_az = _mm512_setzero_pd();

        for (int j = 0; j < N; j += 8) {
            // Cargar 8 elementos
            __m512d xj = _mm512_load_pd(&stars->Cx[j]);
            __m512d yj = _mm512_load_pd(&stars->Cy[j]);
            __m512d zj = _mm512_load_pd(&stars->Cz[j]);
            __m512d mj = _mm512_load_pd(&stars->mass[j]);

            __m512d dx = _mm512_sub_pd(xi, xj);
            __m512d dy = _mm512_sub_pd(yi, yj);
            __m512d dz = _mm512_sub_pd(zi, zj);

            __m512d dist2 = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
            __m512d dist = _mm512_sqrt_pd(dist2);
            __m512d dist3 = _mm512_mul_pd(dist2, dist);

            // Fuerza = -G * m / r^3
            __m512d inv = _mm512_div_pd(G_vec, dist3);
            __m512d F = _mm512_mul_pd(inv, mj);

            // Crear máscara para i != j
            __mmask8 mask = _mm512_cmpneq_epi64_mask(
                _mm512_set1_epi64(i),
                _mm512_set_epi64(j+7, j+6, j+5, j+4, j+3, j+2, j+1, j+0)
            );

            // sum += F * dX con máscara
            sum_ax = _mm512_mask3_fmadd_pd(F, dx, sum_ax, mask);
            sum_ay = _mm512_mask3_fmadd_pd(F, dy, sum_ay, mask);
            sum_az = _mm512_mask3_fmadd_pd(F, dz, sum_az, mask);
        }

        // Reducción horizontal para obtener el valor final
        ax[i] = _mm512_reduce_add_pd(sum_ax);
        ay[i] = _mm512_reduce_add_pd(sum_ay);
        az[i] = _mm512_reduce_add_pd(sum_az);
    }
}
#endif

// Función principal de simulación
void simulate(Star *estrellas, const int N) {
    struct timeval start, end;
    double DT2 = 0.5 * DT;
    double *ax = malloc(N * sizeof(double));
    double *ay = malloc(N * sizeof(double));
    double *az = malloc(N * sizeof(double));
    gettimeofday(&start, NULL);
    for (int step = 0; step < STEPS; step++) {
        compute_aceleration(estrellas, ax, ay, az, N);
        for (int i = 0; i < N; i++) {
            // Leapfrog integration: actualizar velocidad a mitad de paso
            estrellas->Vx[i] = fma(DT2, ax[i], estrellas->Vx[i]);
            estrellas->Vy[i] = fma(DT2, ay[i], estrellas->Vy[i]);
            estrellas->Vz[i] = fma(DT2, az[i], estrellas->Vz[i]);
            // Actualizar posición
            estrellas->Cx[i] = fma(DT, estrellas->Vx[i], estrellas->Cx[i]);
            estrellas->Cy[i] = fma(DT, estrellas->Vy[i], estrellas->Cy[i]);
            estrellas->Cz[i] = fma(DT, estrellas->Vz[i], estrellas->Cz[i]);
        }
        compute_aceleration(estrellas, ax, ay, az, N);
        for (int i = 0; i < N; i++) {
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
    printf("Simuladas %d estrellas en %.2f segundos\n", N, seconds);
    fflush(stdout);
    free(ax);
    free(ay);
    free(az);
    FILE *file = fopen("galactic_positions.txt", "w");
    if (!file) {
        fprintf(stderr, "Error al abrir el archivo\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(file, "ID: %lu X: %.20f Y = %.20f, Z = %.20f\n", estrellas->id[i], estrellas->Cx[i], estrellas->Cy[i],
                estrellas->Cz[i]);
    }
}

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

long octree_new_node(Octree *tree, double cx, double cy, double cz, double half_size) {
    if (tree->size >= tree->capacity) {
        tree->capacity *= 1.5;
        resize_tree(tree);
    }
    size_t i = tree->size++;

    tree->center_x[i] = cx;
    tree->center_y[i] = cy;
    tree->center_z[i] = cz;
    tree->half_size[i] = half_size;

    tree->mass[i] = 0.0;
    tree->com_x[i] = 0.0;
    tree->com_y[i] = 0.0;
    tree->com_z[i] = 0.0;
    tree->star_index[i] = -1;

    for (int j = 0; j < 8; j++)
        tree->children[i][j] = -1;

    return i;
}

inline int get_octant(double cx, double cy, double cz, double x, double y, double z) {
    return ((x >= cx) << 2) | ((y >= cy) << 1) | (z >= cz);
}

void octree_insert(Octree *tree, Star *stars, long node_index, long star_index) {
    double cx = tree->center_x[node_index];
    double cy = tree->center_y[node_index];
    double cz = tree->center_z[node_index];
    double hs = tree->half_size[node_index];

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

    if (tree->children[node_index][oct] == -1) {
        // Crear nuevo nodo hijo
        double offset = hs * 0.5;
        double new_cx = cx + ((oct & 4) ? offset : -offset);
        double new_cy = cy + ((oct & 2) ? offset : -offset);
        double new_cz = cz + ((oct & 1) ? offset : -offset);

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
            if (child != -1) {
                compute_acceleration_bh(stars, tree, child, star_idx, theta, ax, ay, az);
            }
        }
    }
}

void aux_time_bh(const Star *stars, const Octree *tree, long node_idx, long index, double theta, double *ax,
                 double *ay, double *az, double *seconds) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    compute_acceleration_bh(stars, tree, node_idx, index, theta, ax, ay, az);
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

void compute_root_bounds(Star *estrellas, double *center_x, double *center_y, double *center_z, double *half_size) {
    // Inicializar límites
    double min_cx = estrellas->Cx[0], max_cx = estrellas->Cx[0];
    double min_cy = estrellas->Cy[0], max_cy = estrellas->Cy[0];
    double min_cz = estrellas->Cz[0], max_cz = estrellas->Cz[0];

    // Calcular límites de posición
    for (unsigned long i = 0; i < estrellas->size; i++) {
        if (estrellas->Cx[i] < min_cx) min_cx = estrellas->Cx[i];
        else if (estrellas->Cx[i] > max_cx) max_cx = estrellas->Cx[i];
        if (estrellas->Cy[i] < min_cy) min_cy = estrellas->Cy[i];
        else if (estrellas->Cy[i] > max_cy) max_cy = estrellas->Cy[i];
        if (estrellas->Cz[i] < min_cz) min_cz = estrellas->Cz[i];
        else if (estrellas->Cz[i] > max_cz) max_cz = estrellas->Cz[i];
    }

    // Calcular centro
    *center_x = 0.5 * (min_cx + max_cx);
    *center_y = 0.5 * (min_cy + max_cy);
    *center_z = 0.5 * (min_cz + max_cz);

    // Calcular rango máximo
    double dx = max_cx - min_cx;
    double dy = max_cy - min_cy;
    double dz = max_cz - min_cz;
    double max_range = fmax(dx, fmax(dy, dz));

    // Usar margen de seguridad (20%) y dividir entre 2
    *half_size = 0.5 * max_range * 1.2;

    //Elegir precisión para subdivisiones
    MIN_NODE_SIZE = max_range * 1e-7;
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
        for (int j = 0; j < 8; j++) tree->children[i][j] = -1;
        tree->star_index[i] = -1;
    }

    double cx, cy, cz;
    double hs;
    compute_root_bounds(stars, &cx, &cy, &cz, &hs);

    long root = octree_new_node(tree, cx, cy, cz, hs);

    for (unsigned long i = 0; i < stars->size; i++) {
        octree_insert(tree, stars, root, i);
        if (i % 1000000 == 0) {
            printf("\r%ld de %ld estrellas insertadas: %ld nodos en el árbol", i, stars->size, tree->size);
            fflush(stdout);
        }
    }
    if (tree->size < tree->capacity) {
        tree->capacity = tree->size;
        resize_tree(tree);
    }
    gettimeofday(&end, NULL);
    size_t memory = tree->capacity * (
                        sizeof(double) * 7 + // center_x, center_y, center_z, half_size, com_x, com_y, com_z
                        sizeof(double) + // mass
                        sizeof(long[8]) + // children (8 longs por nodo)
                        sizeof(long) // star_index
                    );
    double secs = get_seconds(start, end);
    printf("\nÁrbol de %ld nodos creado en %.4f segundos usando %lu MB \n", tree->capacity, secs, memory / 1024 / 1024);
    fflush(stdout);
    return tree;
}

void test_simulation(Star *estrellas) {
    double THETA[10]= {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    int indexes[20];
    double seconds[20], seconds_bh[200];
    for (int i = 0; i < 20; i++) {
        indexes[i] = rand() % estrellas->size;
    }
    double ax[20] = {0}, ay[20] = {0}, az[20] = {0};
    double axb[200] = {0}, ayb[200] = {0}, azb[200] = {0};

    //free_aux(estrellas);
    Octree *octree = build_tree(estrellas);

    for (int i = 0; i < 20; i++) {
        compute_aceleration_single(estrellas, &ax[i], &ay[i], &az[i], indexes[i], &seconds[i]);
        for (int j = 0; j < 10; j++) {
            int idex=i*10+j;
            aux_time_bh(estrellas, octree, 0, indexes[i], THETA[j], &axb[idex], &ayb[idex], &azb[idex], &seconds_bh[idex]);
        }
        printf("------------------------------------------------------\n");
        printf("Estrella: %d\n", indexes[i]);
        printf("Referencia:              X= %+e Y= %+e Z= %+e  %f segundos\n", ax[i], ay[i], az[i], seconds[i]);
        for (int j = 0; j < 10; j++) {
            int idex=i*10+j;
            printf("BarnesHut THETA %0.1f:     X= %+e Y= %+e Z= %+e  %f segundos\n",THETA[j], axb[idex], ayb[idex], azb[idex],seconds_bh[idex]);
        }
    }

    free_tree(octree);
}