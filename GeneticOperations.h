#ifndef GeneticOperations_H
#define GeneticOperations_H

#include <iostream>
#include <vector>
#include "hamerly/dataset.h"
#include "hamerly/kmeans.h"
#include "hamerly/general_functions.h"
#include "Solution.h"
#include "MathUtils.h"
#include "Param.h"
#include "Hash.h"
#include "Heap.h"
#include "dlib-master/dlib/optimization/max_cost_assignment.h"

using namespace std;

class GeneticOperations {
    
    private:

        
    public:
        
        GeneticOperations();

        Solution* SelectParent(vector<Solution*> pop, int W);

        unsigned short* CreateKmeansSolution(const Dataset* x, unsigned short m);

        unsigned short* CreateKppSolution(const Dataset* x, unsigned short m);

        vector<Solution*> GetInitialPopulation(const Dataset* x, unsigned short m, Param prm);

        int* GetCardinality(int** clusterSize, int m);
        
        vector<Solution*> SelectSurvivors(vector<Solution*> population, const int sizePopulation, Dataset const *x, int m);

        vector<long> MinAssignment(double** c1, double** c2, int m, int d);

        Solution* Crossover(Solution* p1, Solution* p2, const Dataset* x, const int m, double alpha);
};

#endif