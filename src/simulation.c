#include "simulation.h"
#include "aux_fun.h"

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
void compute_aceleration_single(const Star *stars, double *ax, double *ay, double *az, const unsigned long index, double *seconds) {
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

//funcion de prueba: calcula la aceleracion de una sola estrella con las CERCANAS
void compute_aceleration_single_near(const Star *stars, double *ax, double *ay, double *az, const unsigned long index, int *count,
                                     double *seconds) {
    struct timeval start, end;
    const double dx0 = stars->Cx[index];
    const double dy0 = stars->Cy[index];
    const double dz0 = stars->Cz[index];
    const double cutoff = 1000; // pc
    const double cutoff_sq = cutoff * cutoff;
    *count = 0;
    gettimeofday(&start, NULL);
    for (unsigned long i = 0; i < stars->size; i++) {
        if (i != index) {
            double dx = stars->Cx[i] - dx0;
            double dy = stars->Cy[i] - dy0;
            double dz = stars->Cz[i] - dz0;
            double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;

            if (dist_sq < cutoff_sq) {
                (*count)++;
                double dist = sqrt(dist_sq);
                double force = -G * stars->mass[i] / (dist_sq * dist);
                *ax = fma(force, dx, *ax);
                *ay = fma(force, dy, *ay);
                *az = fma(force, dz, *az);
            }
        }
    }
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

void insert_star_to_child(Octree *tree, const Star *stars, long node_idx, long star_idx);

unsigned long octree_add_node(Octree *tree, double cx, double cy, double cz, float half_size) {
    if (tree->size >= tree->capacity) {
        size_t new_capacity = (tree->capacity < 1000000000)
            ? tree->capacity * 2
            : tree->capacity * 1.5;
        tree->capacity = new_capacity;
        resize_tree(tree);
        printf("Nuevo tamaño reservado: %ld\n", tree->capacity);
        fflush(stdout);
    }

    size_t idx = tree->size++;
    tree->cx[idx] = cx;
    tree->cy[idx] = cy;
    tree->cz[idx] = cz;
    tree->half_size[idx] = half_size;
    tree->mass[idx] = 0.0f;
    tree->com_x[idx] = 0.0F;
    tree->com_y[idx] = 0.0F;
    tree->com_z[idx] = 0.0F;
    tree->star_index[idx] = -1;
    tree->first_child_index[idx] = -1;

    return idx;
}


void insert_star(Octree *tree, const Star *stars, unsigned long node_idx, long star_idx) {
    const double x = stars->Cx[star_idx];
    const double y = stars->Cy[star_idx];
    const double z = stars->Cz[star_idx];
    float m_new = stars->mass[star_idx];

    // Si el nodo está vacío, insertar directamente
    if (tree->star_index[node_idx] == -1 && tree->first_child_index[node_idx] == -1) {
        tree->star_index[node_idx] = star_idx;
        tree->mass[node_idx] = m_new;
        tree->com_x[node_idx] = x;
        tree->com_y[node_idx] = y;
        tree->com_z[node_idx] = z;
        return;
    }

    // Si es una hoja con una estrella, subdividir
    if (tree->first_child_index[node_idx] == -1) {
        long existing_star = tree->star_index[node_idx];
        tree->star_index[node_idx] = -1;  // ya no es hoja

        // Crear hijos
        for (int i = 0; i < 8; i++) {
            float offset = tree->half_size[node_idx] / 2.0F;
            double cx = tree->cx[node_idx] + ((i & 1) ? offset : -offset);
            double cy = tree->cy[node_idx] + ((i & 2) ? offset : -offset);
            double cz = tree->cz[node_idx] + ((i & 4) ? offset : -offset);
            unsigned long child_idx = octree_add_node(tree, cx, cy, cz, offset);
            if (i == 0) {
                tree->first_child_index[node_idx] = child_idx;
            }
        }

        // Insertar la estrella previa en el hijo correspondiente
        const double ex = stars->Cx[existing_star];
        const double ey = stars->Cy[existing_star];
        const double ez = stars->Cz[existing_star];
        int eoct = ((ex >= tree->cx[node_idx]) << 0) |
                   ((ey >= tree->cy[node_idx]) << 1) |
                   ((ez >= tree->cz[node_idx]) << 2);
        insert_star(tree, stars, tree->first_child_index[node_idx] + eoct, existing_star);
    }

    // Insertar la nueva estrella en el hijo correspondiente
    int octant = ((x >= tree->cx[node_idx]) << 0) |
                 ((y >= tree->cy[node_idx]) << 1) |
                 ((z >= tree->cz[node_idx]) << 2);
    long child_idx = tree->first_child_index[node_idx] + octant;
    insert_star(tree, stars, child_idx, star_idx);

    // Actualizar centro de masa y masa total
    float m_old = tree->mass[node_idx];
    tree->mass[node_idx] += m_new;
    tree->com_x[node_idx] = (tree->com_x[node_idx] * m_old + x * m_new) / tree->mass[node_idx];
    tree->com_y[node_idx] = (tree->com_y[node_idx] * m_old + y * m_new) / tree->mass[node_idx];
    tree->com_z[node_idx] = (tree->com_z[node_idx] * m_old + z * m_new) / tree->mass[node_idx];
}


void compute_acceleration_bh(const Star *stars, const Octree *tree, unsigned long node_idx, long star_idx,
                             double theta, double *ax, double *ay, double *az) {
    // Evitar autointeracción si es hoja con la misma estrella
    if (tree->first_child_index[node_idx] == -1 && tree->star_index[node_idx] == star_idx) {
        return;
    }

    double dx = tree->com_x[node_idx] - stars->Cx[star_idx];
    double dy = tree->com_y[node_idx] - stars->Cy[star_idx];
    double dz = tree->com_z[node_idx] - stars->Cz[star_idx];
    double dist_sq = dx * dx + dy * dy + dz * dz + EPSILON;
    double dist = sqrt(dist_sq);

    float s = 2.0f * tree->half_size[node_idx];

    // Criterio de Barnes-Hut o si es hoja
    if ((s / dist < theta) || tree->first_child_index[node_idx] == -1) {
        double force = -G * tree->mass[node_idx] / (dist_sq * dist);
        *ax = fma(force, dx, *ax);
        *ay = fma(force, dy, *ay);
        *az = fma(force, dz, *az);
    } else {
        unsigned long base = tree->first_child_index[node_idx];
        for (int i = 0; i < 8; i++) {
            compute_acceleration_bh(stars, tree, base + i, star_idx, theta, ax, ay, az);
        }
    }
}


void aux_time_bh(const Star *stars, const Octree *tree, unsigned long node_idx, long index, double theta, double *ax,
                 double *ay, double *az, double *seconds) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    compute_acceleration_bh(stars, tree, node_idx, index, theta, ax, ay, az);
    gettimeofday(&end, NULL);
    *seconds = get_seconds(start, end);
}

void compute_root_bounds(Star *estrellas, double *center_x, double *center_y, double *center_z, float *half_size) {
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
    *half_size = 0.5F * max_range * 1.2F;
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

    double cx, cy, cz;
    float hs;
    compute_root_bounds(stars, &cx, &cy, &cz, &hs);

    unsigned long root = octree_add_node(tree, cx, cy, cz, hs);
    for (unsigned long i = 0; i < stars->size; i++) {
        insert_star(tree, stars, root, i);
    }
    if (tree->size < tree->capacity) {
        resize_tree(tree);
    }
    gettimeofday(&end, NULL);
    size_t memory = tree->capacity * (sizeof(double) * 3 + // cx, cy, cz
                             sizeof(float) * 5 + // half_size, mass, com_x, com_y, com_z
                             sizeof(long) + // star_index
                             sizeof(long)); // first_child_index
    double secs = get_seconds(start, end);
    printf("Árbol de %ld nodos creado en %.4f segundos usando %lu MB \n", tree->capacity, secs, memory/1024/1024);
    fflush(stdout);
    return tree;

}

void test_simulation(Star *estrellas) {
    int indexes[20];
    int count[20];
    double seconds[20], seconds_near[20], seconds_bh[20];
    for (int i = 0; i < 20; i++) {
        indexes[i] = rand() % estrellas->size;
    }
    double ax[20] = {0}, ay[20] = {0}, az[20] = {0};
    double axn[20] = {0}, ayn[20] = {0}, azn[20] = {0};
    double axb[20] = {0}, ayb[20] = {0}, azb[20] = {0};

    //free_aux(estrellas);
    Octree *octree = build_tree(estrellas);

    for (int i = 0; i < 20; i++) {
        compute_aceleration_single(estrellas, &ax[i], &ay[i], &az[i], indexes[i], &seconds[i]);
        compute_aceleration_single_near(estrellas, &axn[i], &ayn[i], &azn[i], indexes[i], &count[i], &seconds_near[i]);
        aux_time_bh(estrellas, octree, 0, indexes[i], 0.5, &axb[i], &ayb[i], &azb[i], &seconds_bh[i]);
        printf("------------------------------------------------------\n");
        printf("Estrella: %d\n", indexes[i]);
        printf("Aceleracion:                 X=%e Y=%e Z=%e\n", ax[i], ay[i], az[i]);
        printf("Aceleracion cercanas:        X=%e Y=%e Z=%e\n", axn[i], ayn[i], azn[i]);
        printf("Aceleracion Barnes-Hut:      X=%e Y=%e Z=%e\n", axb[i], ayb[i], azb[i]);
        printf("Tiempo:                 %f segundos\n", seconds[i]);
        printf("Tiempo cercanas:        %f segundos\n", seconds_near[i]);
        printf("Tiempo Barnes-Hut:      %f segundos\n", seconds_bh[i]);
    }
    free_tree(octree);
}
