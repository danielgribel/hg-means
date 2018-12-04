#ifndef Evaluator_H
#define Evaluator_H

/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 *
 * Class dedicated to compute external clustering measures, such as
 * Rand, C-Rand, Normalized Mutal Information, and Centroid Index
 */

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <unistd.h>
#include <sys/resource.h>
#include <math.h>
#include "MathUtils.h"
#include "PbData.h"
#include "Solution.h"

using namespace std;

class Evaluator {

    private:

        // Problem Data
        PbData pb_data;

        // Clustering solution produced by the algorithm
        Solution* solution;

        // Clustering solution representing the ground truth
        Solution* ground_truth;

        // Count coefficients for Rand and C-Rand indicators
        void CountRandCoefficients(int& a, int& b, int& c, int& d);
        
    public:

        Evaluator(PbData pb_data, Solution* solution, Solution* ground_truth);
    
        ~Evaluator();

        // Calculate the Rand indicator
        double Rand();

        // Calculate the C-Rand indicator
        double CRand();

        // Calculate the Normalized mutual information indicator
        double Nmi();

        // Calculate the Centroid Index indicator
        double CentroidIndex();
};

#endif