#ifndef Solution_H
#define Solution_H

#include <iostream>
#include <vector>
#include "hamerly/dataset.h"
#include "hamerly/hamerly_kmeans.h"
#include "MathUtils.h"
#include "PbData.h"

#define MUTATION_RATE 0.2

using namespace std;

class Solution {
    
    private:

        // Assignment representation of a solution
        unsigned short* assignment;

        // Centroids representation of a solution
        double** centroids;
        
        // Solution cost, in terms of MSSC objective
        double cost;

        // Alpha weigth
        double alpha;

        // Problem Data
        PbData pb_data;

        // Calculate the centroids representation of a solution (assignment representation supposed to exist)
        void AssignmentToCentroids();

        // Calculate the assignment representation of a solution (centroids representation supposed to exist)
        void CentroidsToAssignment();

        // Initialize assignment data structure
        void InitAssignment();

        // Initialize centroids data structure
        void InitCentroids();

        // Assignment free memory
        void DeleteAssignment() { delete [] assignment; };

        // Centroids free memory
        void DeleteCentroids();

    public:
        
        // Constructor in which we know solution assignment and cost 
        Solution(unsigned short* assignment, double cost, double alpha, PbData pb_data);

        // Constructor in which we know solution centroids and cost 
        Solution(double** centroids, double cost, double alpha, PbData pb_data);

        // Constructor in which we know solution assignment 
        Solution(unsigned short* assignment, double alpha, PbData pb_data);

        // Constructor in which we know solution centroids
        Solution(double** centroids, double alpha, PbData pb_data);

        ~Solution();

        unsigned short* GetAssignment() { return assignment; };

        double** GetCentroids() { return centroids; };

        double GetCost() { return cost; }

        double GetAlpha() { return alpha; }

        // Mutate the solution
        void Mutate();

        // Apply local search to solution (K-means)
        void DoLocalSearch(Dataset const *x);

        // Mutate the alpha weigth
        void MutateAlpha();

        // Repair the solution if it is degenerated (less than m clusters)
        void Repair();
};

#endif