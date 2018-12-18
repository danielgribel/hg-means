
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
#include <math.h>
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

        double crand;

        double nmi;

        double centroid_index;

        // Generate the centroids representation of a solution (assignment representation supposed to exist)
        void AssignmentToCentroids();

        // Generate the assignment representation of a solution (centroids representation supposed to exist)
        void CentroidsToAssignment();

        // Initialize assignment data structure
        void InitAssignment();

        // Initialize centroids data structure
        void InitCentroids();

        // Count coefficients for Rand and C-Rand indicators
        void CountRandCoefficients(Solution* ground_truth, int& a, int& b, int& c, int& d);

        // Remove random center and re-assign data points to closest remaining center
        // Partial assignment with m-1 centers is generated
        void RemoveCenter(int barycenter);

        // Reinsert removed center c in the position of a data point p; re-assign data points to closest center
        // Complete assignment with m centers is generated
        void ReinsertCenter(int c, int p, vector<double> & dist_centroid);

    public:
        
        // Constructor in which solution assignment and cost are known 
        Solution(unsigned short* assignment, double cost, double alpha, PbData pb_data);

        // Constructor in which solution centroids and cost are known
        Solution(vector< vector<double> > & centroids, double cost, double alpha, PbData pb_data);

        // Constructor in which solution assignment is known
        Solution(unsigned short* assignment, double alpha, PbData pb_data);

        // Constructor in which solution centroids are known
        Solution(vector< vector<double> > & centroids, double alpha, PbData pb_data);

        // Assignment free memory
        ~Solution();

        // Get the assignment representation
        unsigned short* GetAssignment() { return assignment; };

        // Get all centroids
        vector< vector<double> > GetCentroids() { return centroids; };

        // Get specific centroid 
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

        // Repair the solution if assignment is degenerated (less than m clusters)
        void Repair();

        /* External Clustering measures */

        // Calculate the C-Rand indicator
        void ComputeCRand(Solution* ground_truth);

        // Calculate the Normalized mutual information indicator
        void ComputeNmi(Solution* ground_truth);

        // Calculate the Centroid Index indicator
        void ComputeCentroidIndex(Solution* ground_truth);

        void ComputeExternalMetrics(unsigned short* assignment);

        double GetCRand() { return crand; };

        double GetNmi() { return nmi; };

        double GetCentroidIndex() { return centroid_index; };
};

#endif