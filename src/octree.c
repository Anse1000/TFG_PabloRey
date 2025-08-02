#include "octree.h"
#include "aux_fun.h"
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "types.h"

double MIN_NODE_SIZE=1e-10;

#ifdef DEBUG_BUILD
// Función para contar estrellas en un árbol recursivamente
long count_stars_in_subtree(Octree *tree, long node_index) {
    if (node_index == INVALID_INDEX) return 0;

    long count = 0;

    // Si este nodo tiene una estrella, contarla
    if (tree->star_index[node_index] >= 0) {
        count = 1;
    }

    // Contar estrellas en los hijos recursivamente
    for (int i = 0; i < 8; i++) {
        count += count_stars_in_subtree(tree, tree->children[node_index][i]);
    }

    return count;
}

// Función para contar estrellas en cada octante del árbol CPU
void count_cpu_octant_stars(Octree *cpu_tree, long *cpu_counts) {
    // Inicializar contadores
    for (int i = 0; i < 8; i++) {
        cpu_counts[i] = 0;
    }

    // Si el árbol está vacío, retornar
    if (cpu_tree->size == 0) return;

    // Contar estrellas en cada octante del nodo raíz del árbol CPU
    for (int octant = 0; octant < 8; octant++) {
        if (cpu_tree->children[0][octant] != INVALID_INDEX) {
            cpu_counts[octant] = count_stars_in_subtree(cpu_tree, cpu_tree->children[0][octant]);
        }
    }
}
#endif

// Función para calcular límites usando centro de masa
void compute_root_bounds_mass_centered(Star *stars, float *cx, float *cy, float *cz, float *hs, double *min_node_size) {
    if (stars->size == 0) return;

    // Calcular centro de masa
    double total_mass = 0.0;
    double com_x = 0.0, com_y = 0.0, com_z = 0.0;

    for (size_t i = 0; i < stars->size; i++) {
        double mass = stars->mass[i];
        total_mass += mass;
        com_x += stars->Cx[i] * mass;
        com_y += stars->Cy[i] * mass;
        com_z += stars->Cz[i] * mass;
    }

    com_x /= total_mass;
    com_y /= total_mass;
    com_z /= total_mass;

    // Encontrar la estrella más lejana del centro de masa
    double max_dist_sq = 0.0;
    for (size_t i = 0; i < stars->size; i++) {
        double dx = stars->Cx[i] - com_x;
        double dy = stars->Cy[i] - com_y;
        double dz = stars->Cz[i] - com_z;
        double dist_sq = dx*dx + dy*dy + dz*dz;
        if (dist_sq > max_dist_sq) {
            max_dist_sq = dist_sq;
        }
    }

    *cx = com_x;
    *cy = com_y;
    *cz = com_z;
    *hs = sqrt(max_dist_sq) * 1.1; // 10% de margen

    *min_node_size = *hs / MIN_SUBDIVISIONS;
    fflush(stdout);
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
    compute_root_bounds_mass_centered(stars, &cx, &cy, &cz, &hs,&MIN_NODE_SIZE);

    long root = octree_new_node(tree, cx, cy, cz, hs);

    for (unsigned long i = 0; i < stars->size; i++) {
        octree_insert(tree, stars, root, i);
    }
    if (tree->size < tree->capacity) {
        tree->capacity = tree->size;
        resize_tree(tree);
    }
    gettimeofday(&end, NULL);

    size_t node_memory = sizeof(double) * 3 + // center_x, center_y, center_z, half_size, com_x, com_y, com_z
                         sizeof(float) + // mass
                         sizeof(float) * 4 +
                         sizeof(unsigned int[8]) + // children (8 longs por nodo)
                         sizeof(long); // star_index
    double secs = get_seconds(start, end);
    printf("Árbol de %ld nodos creado en %.4f segundos usando %lu MB \n", tree->capacity, secs,tree->capacity * node_memory / 1024 / 1024);
#ifdef DEBUG_BUILD
    long cpu_counts[8];
    count_cpu_octant_stars(tree, cpu_counts);
    for (int i=0;i<8;i++) {
        printf("Subarbol %d: %ld nodos -> %lu MB\n",i,cpu_counts[i],cpu_counts[i]*node_memory/1024/1024);
    }
#endif
    fflush(stdout);
    return tree;
}
#ifdef CUDA
Octree **build_tree_gpu(Star *stars,float *center_x,float *center_y, float *center_z) {
    struct timeval start, end;
    size_t initial_capacity = 10000;
    gettimeofday(&start, NULL);
    printf("Iniciando construccion de los subarboles para GPU\n");fflush(stdout);
    
    Octree **trees = malloc(sizeof(Octree*)*8);
    
    // Calcular los límites del espacio total
    float cx, cy, cz, hs;
    compute_root_bounds_mass_centered(stars, &cx, &cy, &cz, &hs, &MIN_NODE_SIZE);
    
    for (int i = 0; i < 8; i++) {
        trees[i] = malloc(sizeof(Octree));
        memset(trees[i], 0, sizeof(Octree));
        trees[i]->capacity = initial_capacity;
        trees[i]->size = 0;
        resize_tree(trees[i]);
        
        // Inicializar todos los nodos como inválidos
        for (size_t j = 0; j < initial_capacity; j++) {
            for (int k = 0; k < 8; k++) trees[i]->children[j][k] = INVALID_INDEX;
            trees[i]->star_index[j] = -1;
        }
        
        // Crear el nodo raíz para este octante
        float offset = hs * 0.5F;
        float oct_cx = cx + ((i & 4) ? offset : -offset);
        float oct_cy = cy + ((i & 2) ? offset : -offset);
        float oct_cz = cz + ((i & 1) ? offset : -offset);
        
        // Crear el nodo raíz con las coordenadas correctas del octante
        long root_index = octree_new_node(trees[i], oct_cx, oct_cy, oct_cz, offset);
        
        // Verificar que el nodo raíz sea efectivamente el índice 0
        if (root_index != 0) {
            printf("Error: el nodo raíz del octante %d no es el índice 0\n", i);
            exit(1);
        }
    }
    
    // Insertar las estrellas en sus respectivos subárboles
    for (unsigned long i = 0; i < stars->size; i++) {
        int octant = get_octant(cx, cy, cz, stars->Cx[i], stars->Cy[i], stars->Cz[i]);
        octree_insert(trees[octant], stars, 0, i);
    }
    
    // Ajustar capacidades
    for (int i = 0; i < 8; i++) {
        if (trees[i]->size < trees[i]->capacity) {
            trees[i]->capacity = trees[i]->size;
            resize_tree(trees[i]);
        }
    }
    
    gettimeofday(&end, NULL);
    
    // Calcular memoria utilizada
    size_t memory[8];
    for (int i = 0; i < 8; i++) {
        memory[i] = trees[i]->capacity * (
                        sizeof(double) * 3 + // center_x, center_y, center_z, half_size, com_x, com_y, com_z
                        sizeof(double) + // mass
                        sizeof(float) * 4 +
                        sizeof(unsigned int[8]) + // children (8 longs por nodo)
                        sizeof(long) // star_index
                    );
    }
    
    printf("Subarboles creados en %.4f ocupando:\n", get_seconds(start, end));
    for (int i = 0; i < 8; i++) {
        printf("Arbol %d: %lu nodos %lu MB\n", i, trees[i]->capacity, memory[i]/1024/1024);
    }
    fflush(stdout);
    
    *center_x = cx;
    *center_y = cy;
    *center_z = cz;
    return trees;
}
#endif