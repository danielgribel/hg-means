#ifndef MathUtils_h
#define MathUtils_h

#include <stdlib.h>
#include <string>
#include <limits>
#include <vector>

namespace MathUtils {
	
    // Infinite float
	const double MAX_FLOAT = std::numeric_limits<double>::max();
	
    // Infinite integer
	const double MAX_INT = std::numeric_limits<int>::max();

    // Given a vector of numbers, find the corresponding index of a given value
    int FindIndex(std::vector<double> values, double key, int first, int last);
    
    // Generate a random number within the given range
    double RandBetween(double min, double max);

    // Probability of a data point to be selected in mutation, given:
    //  -- the current distance to its centroid
    //  -- the sum of distances from all points to their centroids
    //  -- weight alpha
    double Pr(double dist, double sumDist, double alpha, int n);

    // Replace portion of a string
    std::string ReplaceString(std::string s, const std::string &toReplace, const std::string &replaceWith);

    // Squared Euclidean distance of two general data points in R^d
    double SquaredEuclidean(double* a, double* b, int d);

    // Squared Euclidean distance of a point (given by index in dataset) to a centroid
    double PointCenterDist(int p, double* center, int d, double* data);

    // Memory free
    void DeleteMatrix(double** matrix, int m);
}

#endif