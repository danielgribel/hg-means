#include "Solution.h"

Solution::Solution(unsigned short* assignment, double alpha, PbData pb_data) {
	this->assignment = assignment;
    this->alpha = alpha;
    this->pb_data = pb_data;
    InitCentroids();
    AssignmentToCentroids();
}

Solution::Solution(unsigned short* assignment, double cost, double alpha, PbData pb_data) {
	this->assignment = assignment;
    this->cost = cost;
    this->alpha = alpha;
    this->pb_data = pb_data;
    InitCentroids();
    AssignmentToCentroids();
}

Solution::Solution(double** centroids, double cost, double alpha, PbData pb_data) {
	this->centroids = centroids;
    this->cost = cost;
    this->alpha = alpha;
    this->pb_data = pb_data;
    InitAssignment();
    CentroidsToAssignment();
}

Solution::Solution(double** centroids, double alpha, PbData pb_data) {
	this->centroids = centroids;
    this->alpha = alpha;
    this->pb_data = pb_data;
    InitAssignment();
    CentroidsToAssignment();
}

Solution::~Solution() {
    DeleteAssignment();
    DeleteCentroids();
}

void Solution::AssignmentToCentroids() { // O(nd)
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();

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

void Solution::CentroidsToAssignment() { // O(nmd)
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();
    double mindist;
    double dist;
    for(int i = 0; i < n; i++) {
        mindist = MathUtils::MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            dist = MathUtils::pcDist(i, centroids[j], d, data);
            if(dist < mindist) {
                mindist = dist;
                assignment[i] = j;
            }
        }
    }
}

void Solution::Repair() {
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();

    vector<bool> populated(m, false); // the (boolean) list of empty clusters 
    vector<int> listEmpty; // the list of empty clusters
    
    // Identify populated and empty clusters
    for(int i = 0; i < n; i++)
        populated[assignment[i]] = true;

    // Construct the list of empty clusters
    for(int j = 0; j < m; j++) {
        if(populated[j] == false)
            listEmpty.push_back(j);
    }

    if(listEmpty.size() > 0) {
        vector<double> distToCentroid(n); // The distance from data points to their centroids
        double sumDist = 0.0; // Sum of all distances from data points to their centroids
        vector<int> sizes(m, 0); // The size (cardinality) of clusters
        DeleteCentroids();
        InitCentroids();
        AssignmentToCentroids();
        vector<double> pr(n);

        for(int i = 0; i < n; i++) {
            distToCentroid[i] = MathUtils::pcDist(i, centroids[assignment[i]], d, data);
            sumDist = sumDist + distToCentroid[i];
            sizes[assignment[i]] = sizes[assignment[i]]+1;
        }

        double z = 0.0;
        for(int i = 0; i < n; i++) {
            pr[i] = z + MathUtils::Pr(distToCentroid[i], sumDist, alpha, n);
            z = pr[i];
        }

        int it = 0;
        int p;

        while(it < listEmpty.size()) {
            double r = MathUtils::fRand(0.0, pr[n-1]); // O(1)
            int p = MathUtils::findIndex(pr, r, 0, n-1) + 1;
            if(sizes[assignment[p]] > 1) {
                sizes[assignment[p]] = sizes[assignment[p]]-1;
                assignment[p] = listEmpty[it];
                it++;
            }
        }
    }
    DeleteCentroids();
	InitCentroids();
	AssignmentToCentroids();
}

void Solution::Mutate() {
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();
    vector<double> distCentroid(n);
    vector<double> pr (n);
    vector<int> barycenterObj;
    double** newcentroids = new double* [m];
    double mindist;
    double dist;
    int i;

    // Randomly select one centroid to remove it from the solution
    int barycenter = rand() % m;
    double** c3 = centroids;

    // Keep the data points assigned to the centroid to be removed
    for(int i = 0; i < n; i++) { // O(n)
        if(assignment[i] == barycenter) {
            barycenterObj.push_back(i);
        }
    }

    // Re-assign data points to the closest remaining centroid
    for(int q = 0; q < barycenterObj.size(); q++) { // O(|c|md)
        i = barycenterObj[q];
        mindist = MathUtils::MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            if(j != barycenter) {
                dist = MathUtils::pcDist(i, c3[j], d, data);
                if(dist < mindist) {
                    mindist = dist;
                    assignment[i] = j;
                }    
            }
        }
    }

    double sumDist = 0.0;
    for(int i = 0; i < n; i++) {
        distCentroid[i] = MathUtils::pcDist(i, c3[assignment[i]], d, data);
        sumDist = sumDist + distCentroid[i];
    }

    // Get the cumulative distance from each data point $i$ to its centroid defined in $mutation$, i.e., the solution with $m-1$ centroids
    double z = 0.0;
    for(int i = 0; i < n; i++) { // O(nd)
        pr[i] = z + MathUtils::Pr(distCentroid[i], sumDist, alpha, n);
        z = pr[i];
    }

    // Wheel roulette random choose of a data point. data points far from their centroids are more likely to be chosen
    double r = MathUtils::fRand(0.0, pr[n-1]); // O(1)
    int p = MathUtils::findIndex(pr, r, 0, n-1) + 1;
    
    for(int i = 0; i < m; i++) {
        newcentroids[i] = new double[d];
    }

    // Re-insert the removed centroid in the position determined by the roulette wheel
    for(int i = 0; i < m; i++) { // O(md)
        if(i == barycenter) {
            for(int j = 0; j < d; j++) {
                newcentroids[barycenter][j] = data[(p*d)+j];
            }
        } else {
            for(int j = 0; j < d; j++) {
                newcentroids[i][j] = c3[i][j];
            }
        }
    }

    // Re-assign data points according to the solution with $m$ centroids
    for(int i = 0; i < n; i++) { // O(nd)
        if(MathUtils::pcDist(i, newcentroids[barycenter], d, data) < distCentroid[i]) {
            assignment[i] = barycenter;
        }
    }
    Repair();
    MathUtils::deleteMatrix(newcentroids, m);
}

void Solution::MutateAlpha() {
    alpha = alpha + MathUtils::fRand(-MUTATION_RATE, MUTATION_RATE);
    if(alpha > 1.0) {
        alpha = 1.0;
    }
    if(alpha < 0.0) {
        alpha = 0.0;
    }
}

void Solution::DoLocalSearch(Dataset const *x) {
    int m = pb_data.GetM();
    Kmeans *algorithm = new HamerlyKmeans();
    int numThreads = 1;

    // Check for missing initialization
    if(assignment == NULL) {
        cerr << "Please initialize centers first" << endl;
        return;
    }
    if(x == NULL) {
        cerr << "Please load a dataset first" << endl;
        return;
    }

    // Time the execution and get the number of iterations
    algorithm->initialize(x, m, assignment, numThreads);
    int numIt = algorithm->run(MathUtils::MAX_INT);
    cost = algorithm->getSSE();
    DeleteCentroids();
    InitCentroids();
    AssignmentToCentroids();
    delete algorithm;
}

void Solution::InitAssignment() {
	assignment = new unsigned short[pb_data.GetN()];
}

void Solution::InitCentroids() {
    int d = pb_data.GetD();
    int m = pb_data.GetM();
	centroids = new double*[m];
    for(int i = 0; i < m; i++) {
        centroids[i] = new double[d];
        for(int j = 0; j < d; j++) {
            centroids[i][j] = 0.0;
        }
    }
}

void Solution::DeleteCentroids() {
	MathUtils::deleteMatrix(centroids, pb_data.GetM());
}