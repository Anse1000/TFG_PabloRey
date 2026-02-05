#ifndef SIMULATION_H
#define SIMULATION_H
#include "../types.h"
#include "octree_cpu.h"

void simulate(Star *estrellas,int steps, long N, const char* outputfile, float DT);
void test_simulation(Star *estrellas);
void run_full_cpu_validation(Star *estrellas);
void test_theta(Star *stars, int num_samples);
void benchmark_tree_construction(Star *estrellas, int iteraciones);
#endif //SIMULATION_H
