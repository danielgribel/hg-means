#ifndef Evaluator_H
#define Evaluator_H

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <unistd.h>
#include <sys/resource.h>
#include <math.h>
#include "MathUtils.h"

using namespace std;

class Evaluator {

    private:

        int n; // number of data points

        int m; // number of clusters

        int d; // number of dimensions

        unsigned short* y_pred; // assignment representation for solution

        unsigned short* y; // label (class) of data points

        double** c_pred; // solution centroids

        double** c; // true partition centroids

        void countRandCoefficients(int& a, int& b, int& c, int& d);
        
        double ccDist(double* a, double* b);

    public:

        Evaluator(int n, int m, int d, unsigned short* y_pred, unsigned short* y, double** c_pred, double** c);
    
        ~Evaluator();

        double rand();

        double cRand();

        double nmi();

        double centroidIndex();
};

#endif