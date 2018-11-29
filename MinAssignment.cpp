#include "MinAssignment.h"

std::vector<long> minAssignment(double** mat, int m) {
    dlib::matrix<double> cost(m,m);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            cost(i,j) = -1.0 * mat[i][j];
        }
    }
    // Max cost assignment of Dlib with negative costs in matrix
    std::vector<long> assignment = max_cost_assignment(cost);

    return assignment;
}