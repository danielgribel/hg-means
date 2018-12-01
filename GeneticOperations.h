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

        vector<Solution*> population;

        int size_population;

        int max_population;

        int w;

        vector<long> MinAssignment(double** c1, double** c2, int m, int d);
        
        unsigned short* GetKmeansAssignment(const Dataset* x, unsigned short m);

        unsigned short* GetKppAssignment(const Dataset* x, unsigned short m);

        int* GetCardinality(int** clusterSize, int m);

    public:
        
        GeneticOperations(int size_population, int max_population, int w);

        ~GeneticOperations();

        vector<Solution*> GetPopulation() { return population; };

        void AddSolution(Solution* solution) { population.push_back(solution); };

        void DeleteSolution(int i) { delete population[i]; };

        double GetCost(int i) { return population[i]->GetCost(); }

        double GetAlpha(int i) { return population[i]->GetAlpha(); }

        unsigned short* GetAssignment(int i) { return population[i]->GetAssignment(); }

        Solution* SelectParent();

        void CreateInitialPopulation(const Dataset* x, unsigned short m, bool is_mutable);

        void SelectSurvivors(Dataset const *x, int m);

        Solution* Crossover(Solution* p1, Solution* p2, const Dataset* x, const int m);
};

#endif