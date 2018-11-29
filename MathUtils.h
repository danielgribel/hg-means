#ifndef MathUtils_h
#define MathUtils_h

#include <stdlib.h>
#include <string>
#include <limits>
#include <vector>

namespace MathUtils {
	
	const double MAX_FLOAT = std::numeric_limits<double>::max();
	
	const double MAX_INT = std::numeric_limits<int>::max();

    int findIndex(std::vector<double> values, double key, int first, int last);
    
    double fRand(double min, double max);

    double Pr(double dist, double sumDist, double alpha, int n);

    std::string replaceString(std::string s, const std::string &toReplace, const std::string &replaceWith);

    double squaredEuclidean(double* a, double* b, int d);

    double pcDist(int p, double* center, int d, double* data);
}

#endif