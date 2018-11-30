#ifndef Solution_H
#define Solution_H

#include <iostream>
#include <vector>
#include "hamerly/dataset.h"
#include "MathUtils.h"

using namespace std;

class Solution {
    
    private:

        unsigned short* assignment;

        double** centroids;
        
        double cost;

        double alpha;

        int m;

        void initCentroids(int m, int d);

        void deleteAssignment() { delete [] assignment; };

        void deleteMatrix(double** matrix, int m);

    public:
        
        Solution();

        Solution(unsigned short* assignment, double cost, double alpha, const Dataset* x, const int m);

        Solution(double** centroids, double cost, double alpha, const Dataset* x, const int m);

        Solution(unsigned short* assignment, double alpha, const Dataset* x, const int m);

        Solution(double** centroids, double alpha, const Dataset* x, const int m);

        ~Solution();

        unsigned short* getAssignment() { return assignment; };

        double** getCentroids() { return centroids; };

        double getCost() { return cost; }

        double getAlpha() { return alpha; }

        void assignmentToCentroids(const Dataset* x, const int m);

        void centroidsToAssignment(const Dataset* x, const int m);

        void fixSolution(const Dataset* x, const int m, double alpha);
};

#endif