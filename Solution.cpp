#include "Solution.h"

Solution::Solution() {
    this->assignment = NULL;
    this->cost = 0.0;
}

Solution::Solution(unsigned short* assignment, double cost, double alpha) {
    this->assignment = assignment;
    this->cost = cost;
    this->alpha = alpha;
}

Solution::~Solution() {
    deleteAssignment();
}