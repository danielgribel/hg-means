#ifndef MathUtils_h
#define MathUtils_h

/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 *
 * Utility class for general functions
 */

#include <stdlib.h>
#include <string>
#include <limits>
#include <vector>

using namespace std;

namespace MathUtils {
	
    // Infinite float
	const double MAX_FLOAT = std::numeric_limits<double>::max();
	
    // Infinite integer
	const double MAX_INT = std::numeric_limits<int>::max();

    // Given a vector of numbers, find the corresponding index of a given value
    int FindIndex(std::vector<double> values, double key, int first, int last);
    
    // Generate a random number within the given range
    double RandBetween(double min, double max);

    // Probability that a data point is selected in mutation, given:
    //  -- the current distance to its centroid
    //  -- the sum of distances from all points to their centroids
    //  -- weight alpha
    double Pr(double dist, double sumDist, double alpha, int n);

    // Replace portion of a string
    std::string ReplaceString(std::string s, const std::string &toReplace, const std::string &replaceWith);

    // Squared Euclidean distance of two general data points in R^d
    double SquaredEuclidean(vector<double> a, vector<double> b, int d);

    // Squared Euclidean distance of a point (given by its index in the dataset) to a centroid
    double PointCenterDist(int p, vector<double> center, int d, double* data);
}

#endif