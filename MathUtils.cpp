#include "MathUtils.h"

using namespace std;

namespace MathUtils {
    
    int FindIndex(vector<double> & values, double key, int first, int last) {
        if(values[first] <= key && values[first + 1] >= key) {
            return first;
        }

        int imid = first + (last - first)/2;

        if(first == last || imid == first) {
            return -1;
        }

        if(values[imid] > key) {
            return FindIndex(values, key, first, imid);
        } else if(values[imid] <= key) {
            return FindIndex(values, key, imid, last);
        }
    }

    double RandBetween(double min, double max) {
        double f = (double)rand() / RAND_MAX;
        return min + f * (max - min);
    }

    double Pr(double dist, double sumDist, double alpha, int n) {
        return (1.0 * alpha * dist / sumDist) + ((1.0 - alpha)/n);
    }

    string ReplaceString(string s, const string &toReplace, const string &replaceWith) {
        return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
    }

    double SquaredEuclidean(const vector<double> & a, const vector<double> & b, int d) { // O(d)
        double dist = 0.0;
        for(int i = 0; i < d; i++) {
            dist = dist + ( (a[i] - b[i]) * (a[i] - b[i]) );
        }
        return dist;
    }

    double PointCenterDist(int p, vector<double> & center, int d, double *data) { // O(d)
        double dist = 0.0;
        for(int i = 0; i < d; i++) {
			double temp = data[p*d + i] - center[i];
            dist += temp*temp;
        }
        return dist;
    }
}
