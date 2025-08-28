#include "octree.h"
#include "aux_fun.h"
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"

#ifdef DEBUG_BUILD
double test_subdivisions = 1e7;
// Contar nodos en un subárbol recursivamente
long count_nodes_in_subtree(Octree *tree, long node_index) {
    if (node_index == INVALID_INDEX) return 0;

    long count = 1; // Contamos este nodo

    // Contar nodos en los hijos recursivamente
    for (int i = 0; i < 8; i++) {
        count += count_nodes_in_subtree(tree, tree->children[node_index][i]);
    }

    return count;
}

// Contar nodos en cada octante del nodo raíz
void count_cpu_octant_nodes(Octree *cpu_tree, long *cpu_counts) {
    // Inicializar contadores
    for (int i = 0; i < 8; i++) {
        cpu_counts[i] = 0;
    }

    // Si el árbol está vacío, retornar
    if (cpu_tree->size == 0) return;

    // Contar nodos en cada octante del nodo raíz
    for (int octant = 0; octant < 8; octant++) {
        long child_index = cpu_tree->children[0][octant];
        if (child_index != INVALID_INDEX) {
            cpu_counts[octant] = count_nodes_in_subtree(cpu_tree, child_index);
        }
    }
}

void test_tree(Star *stars) {
    Octree *tree = build_tree(stars);
    free_tree(tree);
    test_subdivisions = 1e8;
    tree=build_tree(stars);
    free_tree(tree);
    test_subdivisions = 1e9;
    tree=build_tree(stars);
    free_tree(tree);
    test_subdivisions = 1e10;
    tree=build_tree(stars);
    free_tree(tree);
    test_subdivisions = 1e11;
    tree=build_tree(stars);
    free_tree(tree);
    test_subdivisions = 1e12;
    tree=build_tree(stars);
    free_tree(tree);
}
#endif

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

void octree_insert(Octree *tree, Star *stars, long node_index, long star_index,const float min_node_size) {
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
    if (hs * 2.0 <= min_node_size) {
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

            octree_insert(tree, stars, child, existing_star,min_node_size);
            octree_insert(tree, stars, child, star_index,min_node_size);
        } else {
            octree_insert(tree, stars, child, star_index,min_node_size);
        }
    }
}

Octree *build_tree(Star *stars) {
    struct timeval start, end;
    size_t initial_capacity = 10000;
    gettimeofday(&start, NULL);
#ifdef DEBUG_BUILD
    printf("Iniciando construccion del Arbol\n");
    fflush(stdout);
#endif
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
    float hs,min_node_size;
    compute_root_bounds(stars, &cx, &cy, &cz, &hs,&min_node_size,MIN_SUBDIVISIONS);

    long root = octree_new_node(tree, cx, cy, cz, hs);

    for (unsigned long i = 0; i < stars->size; i++) {
        octree_insert(tree, stars, root, i,min_node_size);
    }
    if (tree->size < tree->capacity) {
        tree->capacity = tree->size;
        resize_tree(tree);
    }
    gettimeofday(&end, NULL);
#ifdef DEBUG_BUILD
    size_t node_memory = sizeof(double) * 3 + // center_x, center_y, center_z, half_size, com_x, com_y, com_z
                     sizeof(float) + // mass
                     sizeof(float) * 4 +
                     sizeof(unsigned int[8]) + // children (8 longs por nodo)
                     sizeof(long); // star_index
    double secs = get_seconds(start, end);
    printf("Árbol de %ld nodos creado en %.4f segundos usando %lu MB\n", tree->capacity, secs,tree->capacity * node_memory / 1024 / 1024);
    long cpu_counts[8];
    count_cpu_octant_nodes(tree, cpu_counts);
    for (int i=0;i<8;i++) {
        printf("Subarbol %d: %ld nodos -> %lu MB\n",i,cpu_counts[i],cpu_counts[i]*node_memory/1024/1024);
    }
    fflush(stdout);
#endif
    return tree;
}

#ifdef CUDA
Octree **build_tree_gpu(Star *stars, const float cx, const float cy, const float cz, const float hs, const float min_node_size,const unsigned int *offsets) {
    struct timeval start, end;
    size_t initial_capacity = 10000;
    gettimeofday(&start, NULL);
#ifdef DEBUG_BUILD
    printf("Iniciando construccion de los subarboles para GPU\n");
    fflush(stdout);
#endif
    Octree **trees = malloc(sizeof(Octree*) * 8);

    // Paralelizar la inicialización de los 8 subárboles
    #pragma omp parallel for num_threads(8)
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
        
        long root_index = octree_new_node(trees[i], oct_cx, oct_cy, oct_cz, offset);
        
        if (root_index != 0) {
            printf("Error: el nodo raíz del octante %d no es el índice 0\n", i);
            exit(1);
        }
    }
    
    // Paralelizar la inserción usando los rangos de reorder_stars
    #pragma omp parallel for num_threads(8)
    for (int octant = 0; octant < 8; octant++) {
        unsigned long start_idx = offsets[octant];
        unsigned long end_idx = (octant == 7) ? stars->size : offsets[octant + 1];
        
        // Insertar todas las estrellas de este octante
        for (unsigned long i = start_idx; i < end_idx; i++) {
            octree_insert(trees[octant], stars, 0, i,min_node_size);
        }
    }
    
    //reajustar tamaños
    for (int i = 0; i < 8; i++) {
        if (trees[i]->size < trees[i]->capacity) {
            trees[i]->capacity = trees[i]->size;
            resize_tree(trees[i]);
        }
    }
    
    gettimeofday(&end, NULL);
#ifdef DEBUG_BUILD
    // Calcular memoria utilizada
    size_t memory[8];
    for (int i = 0; i < 8; i++) {
        memory[i] = trees[i]->capacity * (
                        sizeof(double) * 3 + 
                        sizeof(double) + 
                        sizeof(float) * 4 +
                        sizeof(unsigned int[8]) + 
                        sizeof(long)
                    );
    }
    
    printf("Subarboles creados en %.4f ocupando:\n", get_seconds(start, end));
    for (int i = 0; i < 8; i++) {
        printf("Arbol %d: %lu nodos %lu MB\n", i, trees[i]->capacity, memory[i]/1024/1024);
    }
    fflush(stdout);
#endif
    return trees;
}
#endif