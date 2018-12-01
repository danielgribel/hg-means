#ifndef Solution_H
#define Solution_H

#include <iostream>
#include <vector>
#include "hamerly/dataset.h"
#include "hamerly/hamerly_kmeans.h"
#include "MathUtils.h"

#define MUTATION_RATE 0.2

using namespace std;

class Solution {
    
    private:

        unsigned short* assignment;

        double** centroids;
        
        double cost;

        double alpha;

        double* data;

        int n;

        int d;

        int m;

        void AssignmentToCentroids();

        void CentroidsToAssignment();

        void InitAssignment();

        void InitCentroids();

        void DeleteAssignment() { delete [] assignment; };

        void DeleteCentroids();

    public:
        
        Solution(unsigned short* assignment, double cost, double alpha, const Dataset* x, const int m);

        Solution(double** centroids, double cost, double alpha, const Dataset* x, const int m);

        Solution(unsigned short* assignment, double alpha, const Dataset* x, const int m);

        Solution(double** centroids, double alpha, const Dataset* x, const int m);

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