#include "Solution.h"

Solution::Solution(unsigned short* assignment, double alpha, const Dataset* x, const int m) {
	this->assignment = assignment;
    this->alpha = alpha;
    this->m = m;
    initCentroids(m, x->d);
    assignmentToCentroids(x, m);
}

Solution::Solution(unsigned short* assignment, double cost, double alpha, const Dataset* x, const int m) {
	this->assignment = assignment;
    this->cost = cost;
    this->alpha = alpha;
    this->m = m;
    initCentroids(m, x->d);
    assignmentToCentroids(x, m);
}

Solution::Solution(double** centroids, double cost, double alpha, const Dataset* x, const int m) {
	this->centroids = centroids;
    this->cost = cost;
    this->alpha = alpha;
    this->m = m;
    initAssignment(x->n);
    centroidsToAssignment(x, m);
}

Solution::Solution(double** centroids, double alpha, const Dataset* x, const int m) {
	this->centroids = centroids;
    this->alpha = alpha;
    this->m = m;
    initAssignment(x->n);
    centroidsToAssignment(x, m);
}

Solution::~Solution() {
    deleteAssignment();
    deleteCentroids(centroids, m);
}

void Solution::assignmentToCentroids(const Dataset* x, const int m) { // O(nd)
    const int n = x->n;
    const int d = x->d;
    double *data = x->data;
    vector<int> sizes(m,0);

    for(int i = 0; i < n; i++) {
        sizes[assignment[i]] = sizes[assignment[i]]+1;
        for(int j = 0; j < d; j++) {
            centroids[assignment[i]][j] = centroids[assignment[i]][j] + data[i*d+j];
        }
    }
    
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < d; j++) {
            if(sizes[i] != 0) {
                centroids[i][j] = (1.0*centroids[i][j])/sizes[i];
            }
        }
    }
}

void Solution::centroidsToAssignment(const Dataset* x, const int m) { // O(nmd)
    double mindist;
    double dist;
    for(int i = 0; i < x->n; i++) {
        mindist = MathUtils::MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            dist = MathUtils::pcDist(i, centroids[j], x->d, x->data);
            if(dist < mindist) {
                mindist = dist;
                assignment[i] = j;
            }
        }
    }
}

void Solution::fixSolution(const Dataset* x, const int m, double alpha) {
    vector<bool> populated(m, false); // the (boolean) list of empty clusters 
    vector<int> listEmpty; // the list of empty clusters
    
    // identify populated/empty clusters
    for(int i = 0; i < x->n; i++)
        populated[assignment[i]] = true;

    // construct the list of empty clusters
    for(int j = 0; j < m; j++) {
        if(populated[j] == false)
            listEmpty.push_back(j);
    }

    if(listEmpty.size() > 0) {
        vector<double> distToCentroid(x->n); // the distance from data points to their centroids
        double sumDist = 0.0; // sum of all distances from data points to their centroids
        vector<int> sizes(m, 0); // the size (cardinality) of clusters
        deleteCentroids(centroids, m);
        initCentroids(m, x->d);
        assignmentToCentroids(x, m);
        vector<double> pr(x->n);

        for(int i = 0; i < x->n; i++) {
            distToCentroid[i] = MathUtils::pcDist(i, centroids[assignment[i]], x->d, x->data);
            sumDist = sumDist + distToCentroid[i];
            sizes[assignment[i]] = sizes[assignment[i]]+1;
        }

        double z = 0.0;
        for(int i = 0; i < x->n; i++) {
            pr[i] = z + MathUtils::Pr(distToCentroid[i], sumDist, alpha, x->n);
            z = pr[i];
        }

        int it = 0;
        int p;

        while(it < listEmpty.size()) {
            double r = MathUtils::fRand(0.0, pr[x->n-1]); // O(1)
            int p = MathUtils::findIndex(pr, r, 0, x->n-1) + 1;
            if(sizes[assignment[p]] > 1) {
                sizes[assignment[p]] = sizes[assignment[p]]-1;
                assignment[p] = listEmpty[it];
                it++;
            }
        }
    }
    deleteCentroids(centroids, m);
	initCentroids(m, x->d);
	assignmentToCentroids(x, m);
}

void Solution::initAssignment(int n) {
	assignment = new unsigned short[n];
}

void Solution::initCentroids(int m, int d) {
	centroids = new double*[m];
    for(int i = 0; i < m; i++) {
        centroids[i] = new double[d];
        for(int j = 0; j < d; j++) {
            centroids[i][j] = 0.0;
        }
    }
}

void Solution::deleteCentroids(double** centroids, int m) {
	MathUtils::deleteMatrix(centroids, m);
}