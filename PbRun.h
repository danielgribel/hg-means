#ifndef PbRun_H
#define PbRun_H

#include <iostream>
#include "Solution.h"

class PbRun {
    
    private:

        Solution* solution;

        double time;

        int lastImprovement;

        double nbIter;

        double avgAlpha;

        bool isValid;
    
    public:
        
        PbRun();

        PbRun(Solution* solution, double time, int lastImprovement, double nbIter, double avgAlpha, bool isValid);

        ~PbRun();

        Solution* getSolution() { return solution; };

        double getTime() { return time; }

        int getLastImprovement() { return lastImprovement; }

        double getNbIter() { return nbIter; }

        double getAvgAlpha() { return avgAlpha; }

        bool getIsValid() { return isValid; }

};

#endif