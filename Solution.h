#ifndef Solution_H
#define Solution_H

#include <iostream>

class Solution {
    
    private:

        unsigned short* assignment;
        
        double cost;

        double alpha;

        void deleteAssignment() { delete [] assignment; };

    public:
        
        Solution();

        Solution(unsigned short* assignment, double cost, double alpha);

        ~Solution();

        unsigned short* getAssignment() { return assignment; };

        double getCost() { return cost; }

        double getAlpha() { return alpha; }
};

#endif