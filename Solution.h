
#ifndef Solution_H
#define Solution_H

/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 *
 * The Solution class represents a clustering solution. It can be expressed by two representations:
 * (i) assignment: specifies for each sample the index of the cluster with which it is associated;
 * (ii) centroids: represents the coordinates of the center of each cluster.
 * A Solution object also stores the cost of the solution (in terms of MSSC objective),
 * the mutation parameter alpha and a reference to problem data
 */

#include <iostream>
#include <vector>
#include "hamerly/dataset.h"
#include "hamerly/hamerly_kmeans.h"
#include "MathUtils.h"
#include "PbData.h"

#define MUTATION_RATE 0.2

using namespace std;
using namespace MathUtils;

class Solution {
    
    private:

        // Assignment representation of a solution
        unsigned short* assignment;

        // Centroids representation of a solution
        vector< vector<double> > centroids;
        
        // Solution cost, in terms of MSSC objective
        double cost;

        // Alpha weigth
        double alpha;

        // Problem Data
        PbData pb_data;

        // Generate the assignment representation of a solution (centroids representation supposed to exist)
        void CentroidsToAssignment();

        // Initialize assignment data structure
        void InitAssignment();

        // Initialize centroids data structure
        void InitCentroids();

    public:
        
        // Constructor in which we know solution assignment and cost 
        Solution(unsigned short* assignment, double cost, double alpha, PbData pb_data);

        // Constructor in which we know solution centroids and cost 
        Solution(vector< vector<double> > centroids, double cost, double alpha, PbData pb_data);

        // Constructor in which we know solution assignment 
        Solution(unsigned short* assignment, double alpha, PbData pb_data);

        // Constructor in which we know solution centroids
        Solution(vector< vector<double> > centroids, double alpha, PbData pb_data);
        
        // Assignment free memory
        ~Solution();

        // Get the assignment representation
        unsigned short* GetAssignment() { return assignment; };

        // Get all centroids
        vector< vector<double> > GetCentroids() { return centroids; };

        // Get especific centroid 
        vector<double> GetCentroids(int i) { return centroids[i]; };

        // Get solution cost (in terms of MSSC objective)
        double GetCost() { return cost; }

        // Get alpha parameter
        double GetAlpha() { return alpha; }

        // Mutate the solution
        void Mutate();

        // Apply local search to solution (K-means)
        void DoLocalSearch(Dataset const *x);

        // Mutate the alpha weigth
        void MutateAlpha();

        // Repair the solution if it is degenerated (less than m clusters)
        void Repair();

        // Generate the centroids representation of a solution (assignment representation supposed to exist)
        void AssignmentToCentroids();
};

#endif