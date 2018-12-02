#ifndef Evaluator_H
#define Evaluator_H

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

        PbData pb_data;

        Solution* solution;

        Solution* ground_truth;

        void CountRandCoefficients(int& a, int& b, int& c, int& d);
        
    public:

        Evaluator(PbData pb_data, Solution* solution, Solution* ground_truth);
    
        ~Evaluator();

        double Rand();

        double CRand();

        double Nmi();

        double CentroidIndex();
};

#endif