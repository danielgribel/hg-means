#ifndef GeneticOperations_H
#define GeneticOperations_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "hamerly/dataset.h"
#include "hamerly/kmeans.h"
#include "hamerly/general_functions.h"
#include "PbData.h"
#include "Solution.h"
#include "MathUtils.h"
#include "Hash.h"
#include "dlib-master/dlib/optimization/max_cost_assignment.h"

using namespace std;

class GeneticOperations {
    
    private:

        vector<Solution*> population;

        Solution* best_solution;

        PbData pb_data;

        Param param;

        vector<long> MinAssignment(double** c1, double** c2);
        
        unsigned short* GetKmeansAssignment(const Dataset* x);

        unsigned short* GetKppAssignment(const Dataset* x);

        int* GetCardinality(int** clusterSize);

        void PushMax(vector< pair<double, int> >& heap, double cost, int val);
        
        int PopMax(vector< pair<double, int> >& heap);
        
        pair<double, int> FrontMax(vector< pair<double, int> >& heap);

    public:
        
        GeneticOperations(PbData pb_data, Param param);

        ~GeneticOperations();

        vector<Solution*> GetPopulation() { return population; };

        void AddSolution(Solution* solution) { population.push_back(solution); };

        void DeleteSolution(int i) { delete population[i]; };

        double GetCost(int i) { return population[i]->GetCost(); }

        double GetAlpha(int i) { return population[i]->GetAlpha(); }

        unsigned short* GetAssignment(int i) { return population[i]->GetAssignment(); }

        Solution* SelectParent();

        void CreateInitialPopulation(const Dataset* x);

        void SelectSurvivors(const Dataset* x);

        Solution* Crossover(Solution* p1, Solution* p2);

        Solution* GetBestSolution() { return best_solution; };

        void StoreBestSolution(Solution* s);

        void ReplaceBestSolution(Solution* s);

        void HGMeans(const Dataset* x);
};

#endif