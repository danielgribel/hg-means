#ifndef PbRun_H
#define PbRun_H

#include <iostream>
#include "Solution.h"

class PbRun {
    
    private:

        Solution* solution;

        double time;
    
    public:
        
        PbRun();

        PbRun(Solution* solution, double time);

        ~PbRun();

        Solution* GetSolution() { return solution; };

        double GetTime() { return time; }
};

#endif