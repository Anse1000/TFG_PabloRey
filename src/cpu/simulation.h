#ifndef SIMULATION_H
#define SIMULATION_H
#include "../types.h"
#include "octree_cpu.h"

void simulate(Star *estrellas,int steps, long N, const char* outputfile, float DT);
void test_simulation(Star *estrellas);
#endif //SIMULATION_H
