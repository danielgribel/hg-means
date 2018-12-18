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

struct Item {
    double cost;
    vector<int> cardinality; // must be a sorted vector
};

struct KeyHash {
    public:
    long operator () ( const Item& x ) const {
        long seed = 0;
        for(int i = 0; i < x.cardinality.size(); i++) {
            seed = seed + (i+1)*x.cardinality[i];
        }
        return long(x.cost + seed) % 99991;
    }
};

struct KeyEqual {
    public:
    bool operator () ( const Item& x, const Item& y ) const {
        if ( !(x.cost > y.cost - 0.00000001 && x.cost < y.cost + 0.00000001) ) {
            return false;
        }
        for(int i = 0; i < x.cardinality.size(); i++) {
            if(x.cardinality[i] != y.cardinality[i]) {
                return false;
            }
        }
		return true ;
    }
};

namespace MathUtils {
	
    // Infinite float
	const double MAX_FLOAT = std::numeric_limits<double>::max();
	
    // Infinite integer
	const int MAX_INT = std::numeric_limits<int>::max();

    // Given a vector of numbers, find the corresponding index of a given value
    int FindIndex(std::vector<double> & values, double key, int first, int last);
    
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
    double SquaredEuclidean(const vector<double> & a, const vector<double> & b, int d);

    // Squared Euclidean distance of a point (given by its index in the dataset) to a centroid
    double PointCenterDist(int p, vector<double> & center, int d, double* data);
}

#endif