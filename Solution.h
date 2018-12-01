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

        unsigned short* assignment;

        double** centroids;
        
        double cost;

        double alpha;

        PbData pb_data;

        void AssignmentToCentroids();

        void CentroidsToAssignment();

        void InitAssignment();

        void InitCentroids();

        void DeleteAssignment() { delete [] assignment; };

        void DeleteCentroids();

    public:
        
        Solution(unsigned short* assignment, double cost, double alpha, PbData pb_data);

        Solution(double** centroids, double cost, double alpha, PbData pb_data);

        Solution(unsigned short* assignment, double alpha, PbData pb_data);

        Solution(double** centroids, double alpha, PbData pb_data);

        ~Solution();

        unsigned short* GetAssignment() { return assignment; };

        double** GetCentroids() { return centroids; };

        double GetCost() { return cost; }

        double GetAlpha() { return alpha; }

        void Mutate();

        void DoLocalSearch(Dataset const *x);

        void MutateAlpha();

        void Repair();
};

#endif