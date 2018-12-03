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

        // Population of solutions
        vector<Solution*> population;

        // The best solution in the population
        Solution* best_solution;

        // Problem data
        PbData pb_data;

        // Problem parameters
        Param param;

        // Get the minimum assignment between two centroids. Hungarian method of Dlib is used
        vector<long> MinAssignment(vector< vector<double> > c1, vector< vector<double> > c2);
        
        // Generate a solution for MSSC with initial centers randomly chosen 
        unsigned short* GetKmeansAssignment(const Dataset* x);
        
        // Generate a solution for MSSC with initial centers chosen from K-means++ heuristic
        unsigned short* GetKppAssignment(const Dataset* x);

        // Get the cardinalities (number of points) of clusters
        int* GetCardinality(int** clusterSize);

        // Push element to heap, such that the element with maximum value is on top
        void PushMax(vector< pair<double, int> >& heap, double cost, int val);
        
        // Pop the maximum-valued element from heap
        int PopMax(vector< pair<double, int> >& heap);
        
        // Get the maximum-valued element from heap
        pair<double, int> FrontMax(vector< pair<double, int> >& heap);

    public:
        
        GeneticOperations(PbData pb_data, Param param);

        ~GeneticOperations();

        PbData GetPbData() { return pb_data; };

        Param GetParam() { return param; };

        // Get the population of solutions
        vector<Solution*> GetPopulation() { return population; };

        // Add solution to current population
        void AddSolution(Solution* solution) { population.push_back(solution); };

        // Delete solution to current population
        void DeleteSolution(int i) { delete population[i]; };

        // Get the MSSC cost of a solution
        double GetCost(int i) { return population[i]->GetCost(); }

        // Get the alpha value of a solution
        double GetAlpha(int i) { return population[i]->GetAlpha(); }

        // Get the assignment of a solution
        unsigned short* GetAssignment(int i) { return population[i]->GetAssignment(); }

        // Parent selection: W-tournament selection among the solutions of the population
        Solution* SelectParent();

        // Generate the initial population
        void CreateInitialPopulation(const Dataset* x);

        // Select the survivor solutions for the next generation
        void SelectSurvivors(const Dataset* x);

        // Perform exact minimum-cost matching crossover
        Solution* Crossover(Solution* p1, Solution* p2);

        // Get the best solution of the population -- in terms of MSSC objective
        Solution* GetBestSolution() { return best_solution; };

        // Store the best solution of the population
        void StoreBestSolution(Solution* s);

        // Delete current best solution and replace it
        void ReplaceBestSolution(Solution* s);

        // Execute the Hybrid Genetic algorithm for the MSSC
        void HGMeans(const Dataset* x);
};

#endif