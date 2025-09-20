#include "octree_gpu.h"
#include <stdlib.h>
#include "../aux_fun.h"
#include "../types.h"

void resize_octant(OctreeOctant *tree) {
    tree->nodes = safe_realloc(tree->nodes, sizeof(OctreeNode) * tree->capacity);
}
void resize_frontier(FrontierGPU *tree) {
    tree->nodes = safe_realloc(tree->nodes, sizeof(FrontierNode) * tree->capacity);
}

void free_tree_gpu(OctreeGPU *tree) {
#ifdef DEBUG_BUILD
    printf("Liberando memoria del arbol\n");
    fflush(stdout);
#endif
    if (!tree) return;
    for (int i = 0; i < 8; i++) {
        if (tree->octants[i]) {
            free(tree->octants[i]->nodes);
            free(tree->octants[i]);
        }
        if (tree->frontier[i]) {
            free(tree->frontier[i]->nodes);
            free(tree->frontier[i]);
        }
    }
    free(tree);
#ifdef DEBUG_BUILD
    printf("Memoria liberada\n");
    fflush(stdout);
#endif
}

unsigned int new_node(OctreeOctant *tree,float cx, float cy, float cz, float hs,float m,double com_x,double com_y, double com_z) {
    if (tree->size >= tree->capacity) {
        tree->capacity *= 1.2;
        resize_octant(tree);
    }
    size_t i = tree->size++;

    tree->nodes[i].center_x = cx;
    tree->nodes[i].center_y = cy;
    tree->nodes[i].center_z = cz;
    tree->nodes[i].half_size = hs;
    tree->nodes[i].mass = m;
    tree->nodes[i].com_x = com_x;
    tree->nodes[i].com_y = com_y;
    tree->nodes[i].com_z = com_z;
    tree->nodes[i].star_index = -1;
    tree->nodes[i].next = INVALID_INDEX;
    for (int j = 0; j < 8; j++)
        tree->nodes[i].children[j] = INVALID_INDEX;
    return i;
}

void octree_insert(OctreeOctant *tree, Star *stars,unsigned int node_index, long star_index,const float min_node_size) {
    OctreeNode *node = &tree->nodes[node_index];
    double x = stars->Cx[star_index];
    double y = stars->Cy[star_index];
    double z = stars->Cz[star_index];
    float m = stars->mass[star_index];
    // Actualizar masa y centro de masa del nodo
    float old_mass = node->mass;
    float new_mass = old_mass + m;
    node->com_x = (node->com_x * old_mass + x * m) / new_mass;
    node->com_y = (node->com_y * old_mass + y * m) / new_mass;
    node->com_z = (node->com_z * old_mass + z * m) / new_mass;
    node->mass = new_mass;

    if (node->half_size * 2.0 <= min_node_size) {
        if (node->star_index == -1) {
            node->star_index = star_index;
        }
        return;
    }
    int oct = get_octant(node->center_x, node->center_y, node->center_z, x, y, z);

    if (node->children[oct] == INVALID_INDEX) {
        float offset = node->half_size * 0.5F;
        float new_cx = node->center_x + ((oct & 4) ? offset : -offset);
        float new_cy = node->center_y + ((oct & 2) ? offset : -offset);
        float new_cz = node->center_z + ((oct & 1) ? offset : -offset);
        unsigned int child_index = new_node(tree, new_cx, new_cy, new_cz, offset,m,x,y,z);
        tree->nodes[node_index].children[oct] = child_index;
        tree->nodes[child_index].star_index = star_index;
    }else {
        unsigned int child_index = node->children[oct];
        OctreeNode *child = &tree->nodes[child_index];
        if (child->star_index >= 0) {
            long existing_star = child->star_index;
            child->star_index = -1;
            child->mass = 0.0F;
            child->com_x = child->com_y = child->com_z = 0.0;

            octree_insert(tree, stars, child_index, existing_star,min_node_size);
            octree_insert(tree, stars, child_index, star_index,min_node_size);
        } else {
            octree_insert(tree, stars, child_index, star_index,min_node_size);
        }
    }

}

//inicializar raiz de cada octante
void init_octant(OctreeOctant *tree, float cx, float cy, float cz, float half_size) {
    tree->size=1;
    OctreeNode *n = &tree->nodes[0];
    n->center_x = cx;
    n->center_y = cy;
    n->center_z = cz;
    n->half_size = half_size;
    n->mass = 0.0f;
    n->com_x = n->com_y = n->com_z = 0.0;
    n->star_index = -1;
    n->next = INVALID_INDEX;
    for (int i=0;i<8;i++) n->children[i] = INVALID_INDEX;
}
void build_frontier_from_octant(const OctreeOctant *tree, unsigned int node_idx,
                                float cx, float cy, float cz,
                                float half_size,float theta2, FrontierGPU *frontier)
{
    const OctreeNode *node = &tree->nodes[node_idx];

    // distancia al centro del octante destino
    double dx = node->com_x - cx;
    double dy = node->com_y - cy;
    double dz = node->com_z - cz;
    double dist2 = dx*dx + dy*dy + dz*dz + 1e-12;

    double s = 2.0 * node->half_size;
    double s2 = s * s;

    // Criterio Barnes–Hut: ¿nodo suficientemente lejano?
    if (s2 < theta2 * dist2 || node->star_index >= 0) {
        // aceptar este nodo como partícula frontera
        FrontierNode fn;
        fn.com_x = node->com_x;
        fn.com_y = node->com_y;
        fn.com_z = node->com_z;
        fn.mass  = node->mass;

        if (frontier->size >= frontier->capacity) {
            frontier->capacity *= 1.2;
            resize_frontier(frontier);
        }
        size_t i = frontier->size++;
        frontier->nodes[i] = fn;
    } else {
        // bajar a los hijos
        for (int i=0;i<8;i++) {
            unsigned int c = node->children[i];
            if (c != INVALID_INDEX) {
                build_frontier_from_octant(tree, c,
                                           cx, cy, cz, half_size,
                                           theta2, frontier);
            }
        }
    }
}

void build_frontier(OctreeGPU *tree,int octant) {
    FrontierGPU *frontier = tree->frontier[octant];
    float hs = tree->octants[octant]->nodes[0].half_size;
    float cx = tree->octants[octant]->nodes[0].center_x;
    float cy = tree->octants[octant]->nodes[0].center_y;
    float cz = tree->octants[octant]->nodes[0].center_z;
    float theta2 = THETA * THETA;

    for (int i = 0; i < 8; i++) {
        if (i==octant) continue;
        build_frontier_from_octant(tree->octants[i],0,cx,cy,cz,hs,theta2,frontier);
    }
    if (frontier->size < frontier->capacity) {
        frontier->capacity = frontier->size;
        resize_frontier(frontier);
    }

}
void build_ropes(OctreeNode *nodes, unsigned int node_idx, unsigned int next) {
    if (node_idx == INVALID_INDEX) return;

    // asignar rope al nodo actual
    nodes[node_idx].next = next;

    // recorrer hijos en orden
    for (int i = 0; i < 8; i++) {
        unsigned int child = nodes[node_idx].children[i];
        if (child == INVALID_INDEX) continue;

        // calcular el rope para el hijo:
        unsigned int child_next = INVALID_INDEX;

        // buscar siguiente hermano válido
        for (int j = i + 1; j < 8; j++) {
            if (nodes[node_idx].children[j] != INVALID_INDEX) {
                child_next = nodes[node_idx].children[j];
                break;
            }
        }

        // si no hay hermano, usar el rope del padre
        if (child_next == INVALID_INDEX) {
            child_next = next;
        }
        // llamada recursiva
        build_ropes(nodes, child, child_next);
    }
}
OctreeGPU *build_tree(Star *stars, const float cx, const float cy, const float cz, const float hs, const float min_node_size,const unsigned int *offsets) {
    struct timeval start, end;
    size_t initial_capacity = 10000;
    gettimeofday(&start, NULL);
#ifdef DEBUG_BUILD
    printf("Iniciando construccion del Arbol\n");
    fflush(stdout);
#endif

    OctreeGPU *tree = malloc(sizeof(OctreeGPU));
    for (int i = 0; i < 8; i++) {
        tree->octants[i] = malloc(sizeof(OctreeOctant));
        tree->frontier[i] = malloc(sizeof(FrontierGPU));
        //inicializar octante
        tree->octants[i]->capacity = initial_capacity;
        tree->octants[i]->size = 0;
        tree->octants[i]->nodes = NULL;
        resize_octant(tree->octants[i]);

        // Inicializar frontera
        tree->frontier[i]->capacity = 1000;
        tree->frontier[i]->size = 0;
        tree->frontier[i]->nodes = NULL;
        resize_frontier(tree->frontier[i]);
    }

 #pragma omp parallel for num_threads(8)
    for (int i = 0; i < 8; i++) {
        OctreeOctant *octant = tree->octants[i];

        float offset = hs * 0.5F;
        float oct_cx = cx + ((i & 4) ? offset : -offset);
        float oct_cy = cy + ((i & 2) ? offset : -offset);
        float oct_cz = cz + ((i & 1) ? offset : -offset);

        init_octant(octant, oct_cx, oct_cy, oct_cz, offset);

        //construir octante
        long start_idx = offsets[i];
        long end_idx = (i == 7) ? stars->size : offsets[i + 1];

        // Insertar todas las estrellas de este octante
        for (long j = start_idx; j < end_idx; j++) {
            octree_insert(octant, stars, 0, j,min_node_size);
        }
        //construir ropes
        build_ropes(octant->nodes,0,INVALID_INDEX);
        //reajustar tamaño
        if (octant->size < octant->capacity) {
            octant->capacity = octant->size;
            resize_octant(octant);
        }
    }
    for (int i = 0; i < 8; i++) {
        //construir frontera
        build_frontier(tree,i);
    }
    gettimeofday(&end, NULL);
#ifdef DEBUG_BUILD
    // Calcular memoria utilizada
    size_t memory[8];
    for (int i = 0; i < 8; i++) {
        memory[i] = tree->octants[i]->capacity * sizeof(OctreeNode);
    }
    printf("Subarboles creados en %.4f segundos ocupando:\n", get_seconds(start, end));
    for (int i = 0; i < 8; i++) {
        printf("Arbol %d: %lu nodos %lu MB\n", i, tree->octants[i]->capacity, memory[i]/1024/1024);
    }
    fflush(stdout);
#endif
    return tree;
}