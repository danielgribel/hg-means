#include "Solution.h"

Solution::Solution(unsigned short* assignment, double alpha, PbData pb_data) {
	this->assignment = assignment;
    this->alpha = alpha;
    this->pb_data = pb_data;
    AssignmentToCentroids();
}

Solution::Solution(unsigned short* assignment, double cost, double alpha, PbData pb_data) {
	this->assignment = assignment;
    this->cost = cost;
    this->alpha = alpha;
    this->pb_data = pb_data;
    AssignmentToCentroids();
}

Solution::Solution(vector< vector<double> > & centroids, double cost, double alpha, PbData pb_data) {
	this->centroids = centroids;
    this->cost = cost;
    this->alpha = alpha;
    this->pb_data = pb_data;
    CentroidsToAssignment();
}

Solution::Solution(vector< vector<double> > & centroids, double alpha, PbData pb_data) {
	this->centroids = centroids;
    this->alpha = alpha;
    this->pb_data = pb_data;
    CentroidsToAssignment();
}

Solution::~Solution() {
    delete [] assignment;
}

void Solution::InitAssignment() {
	assignment = new unsigned short[pb_data.GetN()];
}

void Solution::InitCentroids() {
    centroids = vector< vector<double> > (pb_data.GetM(), vector<double>(pb_data.GetD(), 0.0));
}

void Solution::AssignmentToCentroids() { // O(nd)
    InitCentroids();
    int d = pb_data.GetD();
    double* data = pb_data.GetData();
    vector<int> sizes(pb_data.GetM(),0);

    for(int i = 0; i < pb_data.GetN(); i++) {
        sizes[ assignment[i] ]++;
        for(int j = 0; j < d; j++) {
            centroids[ assignment[i] ][j] = centroids[ assignment[i] ][j] + data[i*d+j];
        }
    }
    for(int i = 0; i < pb_data.GetM(); i++) {
        for(int j = 0; j < d; j++) {
            if(sizes[i] != 0) {
                centroids[i][j] = (1.0*centroids[i][j])/sizes[i];
            }
        }
    }
}

void Solution::CentroidsToAssignment() { // O(nmd)
    InitAssignment();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();
    double dist, mindist;

    for(int i = 0; i < pb_data.GetN(); i++) {
        mindist = MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            dist = PointCenterDist(i, centroids[j], d, data);
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

    // Size (cardinality) of clusters
    vector<int> sizes(m, 0);

    // List of empty clusters (index of clusters)
    vector<int> empty_clusters;
    
    // Count cardinality of each cluster
    for(int i = 0; i < n; i++) {
        sizes[ assignment[i] ]++;
    }

    // Construct the list of empty clusters
    for(int j = 0; j < m; j++) {
        if(sizes[j] == 0)
            empty_clusters.push_back(j);
    }

    if(empty_clusters.size() > 0) {
        vector<double> dist_centroid(n); // Distance from data points to their centroids
        double total_dist = 0.0; // Total distances from data points to their centroids
        vector<double> pr(n);

        // Update the coordinates of centroids
        AssignmentToCentroids();

        // Get the distance from each data point to its centroid
        for(int i = 0; i < n; i++) {
            dist_centroid[i] = PointCenterDist(i, centroids[ assignment[i] ], d, data);
            total_dist = total_dist + dist_centroid[i];
        }

        // Get for each data point the cumulative probability of being chosen to move to an empty center
        pr[0] = Pr(dist_centroid[0], total_dist, alpha, n);
        for(int i = 1; i < n; i++) {
            pr[i] = pr[i-1] + Pr(dist_centroid[i], total_dist, alpha, n);
        }

        int e = 0;
        // While an empty cluster exists, randomly select a data point and move it to an empty cluster
        while(e < empty_clusters.size()) {
            // Wheel roulette selection
            // Portion associated to data points far from their centroids are more likely to be moved to an empty cluster, according to alpha
            double r = RandBetween(0.0, pr[n-1]); // O(1)
            // Find in the wheel roulette the corresponding data point index
            int p = FindIndex(pr, r, 0, n-1) + 1;
            // If the current cluster of point p has more than 1 element, move p to an empty cluster  
            if(sizes[ assignment[p] ] > 1) {
                // Decrement the size of the current cluster of point p
                sizes[ assignment[p] ]--;
                // Move point p to the e-th empty cluster
                assignment[p] = empty_clusters[e];
                e++;
            }
        }
    }
    // Update the coordinates of centroids
	AssignmentToCentroids();
}

void Solution::RemoveCenter(int c) {
    double* data = pb_data.GetData();
    vector<int> list_points;
    double dist, mindist;
    
    // Keep the list of data points assigned to the centroid to be removed
    for(int i = 0; i < pb_data.GetN(); i++) {
        if(assignment[i] == c) {
            list_points.push_back(i);
        }
    }

    // Re-assign data points to the closest remaining centroid
    for(int i = 0; i < list_points.size(); i++) { // O( |c| m d )
        mindist = MAX_FLOAT;
        for(int j = 0; j < pb_data.GetM(); j++) {
            if(j != c) {
                dist = PointCenterDist(list_points[i], centroids[j], pb_data.GetD(), data);
                if(dist < mindist) {
                    mindist = dist;
                    assignment[ list_points[i] ] = j;
                }
            }
        }
    }
}

void Solution::ReinsertCenter(int c, int p, vector<double> & dist_centroid) {
    int d = pb_data.GetD();
    double* data = pb_data.GetData();
    
    // Place centroid c in the position of data point p
    for(int j = 0; j < d; j++) {
        centroids[c][j] = data[(p*d)+j];
    }

    // Re-assign data points according to the solution with m centroids
    for(int i = 0; i < pb_data.GetN(); i++) { // O(nd)
        if(PointCenterDist(i, centroids[c], d, data) < dist_centroid[i]) {
            assignment[i] = c;
        }
    }
}

void Solution::Mutate() {
    int n = pb_data.GetN();
    double* data = pb_data.GetData();
    vector<double> dist_centroid (n);
    vector<double> pr (n);
    double total_dist = 0.0;

    // Randomly select one centroid to remove from the solution
    int c = rand() % pb_data.GetM();

    // Remove random center and re-assign data points to the closest remaining center
    RemoveCenter(c);

    // Get the distance from each data point to its centroid defined in the solution with m-1 centroids
    for(int i = 0; i < n; i++) {
        dist_centroid[i] = PointCenterDist(i, centroids[ assignment[i] ], pb_data.GetD(), data);
        total_dist = total_dist + dist_centroid[i];
    }

    // Get for each data point the cumulative probability of being chosen
    pr[0] = Pr(dist_centroid[0], total_dist, alpha, n);
    for(int i = 1; i < n; i++) {
        pr[i] = pr[i-1] + Pr(dist_centroid[i], total_dist, alpha, n);
    }

    // Wheel roulette selection
    // Portion associated to data points far from their centroids are more likely to be chosen, according to alpha
    double r = RandBetween(0.0, pr[n-1]); // O(1)

    // Find in the wheel roulette the corresponding data point index
    int p = FindIndex(pr, r, 0, n-1) + 1;

    // Reinsert removed center c in the position of a data point p; re-assign data points to closest center
    ReinsertCenter(c, p, dist_centroid);

    // Repair the solution if assignment is degenerated
    Repair();
}

void Solution::MutateAlpha() {
    // Perturbate alpha parameter
    alpha = alpha + RandBetween(-MUTATION_RATE, MUTATION_RATE);
    
    // Force alpha to be within the interval 
    if(alpha > 1.0) {
        alpha = 1.0;
    }
    if(alpha < 0.0) {
        alpha = 0.0;
    }
}

void Solution::DoLocalSearch(Dataset const *x) {
    Kmeans *algorithm = new HamerlyKmeans();

    // Check for missing initialization
    if(assignment == NULL) {
        cerr << "Please initialize centers first" << endl;
        return;
    }
    if(x == NULL) {
        cerr << "Please load a dataset first" << endl;
        return;
    }
    // Initialize centers
    algorithm->initialize(x, pb_data.GetM(), assignment, 1);
    
    // Run K-means
    int numIt = algorithm->run(MAX_INT);

    // Get solution cost
    cost = algorithm->getSSE();
    
    // Update centers
    AssignmentToCentroids();
    
    delete algorithm;
}

void Solution::CountRandCoefficients(Solution* ground_truth, int& a, int& b, int& c, int& d) {
	a = 0;
	b = 0;
	c = 0;
	d = 0;

	unsigned short* y = ground_truth->GetAssignment();

	for(int i = 0; i < pb_data.GetN(); i++) {
		for(int j = i+1; j < pb_data.GetN(); j++) {
			if( (assignment[i] == assignment[j]) && (y[i] == y[j]) ) {
				a++;
			}
			if( (assignment[i] == assignment[j]) && (y[i] != y[j]) ) {
				b++;
			}
			if( (assignment[i] != assignment[j]) && (y[i] == y[j]) ) {
				c++;
			}
			if( (assignment[i] != assignment[j]) && (y[i] != y[j]) ) {
				d++;
			}
		}
	}
}

void Solution::ComputeCRand(Solution* ground_truth) {
	int a, b, c, d;
	CountRandCoefficients(ground_truth, a, b, c, d);
	int total = a + b + c + d;
	double crandIndex = (a - (1.0*(b + a)*(c + a))/total)/((1.0*(b + a + c + a))/2 - (1.0*(b + a)*(c + a))/total);
	crand = crandIndex;
}

// Implemented by @Carlo Nicolini
// More information in the original repository: https://github.com/CarloNicolini/rnmi
void Solution::ComputeNmi(Solution* ground_truth) {
	int n = pb_data.GetN();
	unsigned short* pa = assignment;
	unsigned short* pb = ground_truth->GetAssignment();
	int qa = -1;
	int qb = -1;
	vector <int > ga; // Group A
	vector <int > gb; // Group B
	for(int i = 0; i < n; i++) {
		if(qa < pa[i]) qa = pa[i];
		if(qb < pb[i]) qb = pb[i];
	}
	qa++;
	qb++;

	if(qa == 1 && qb == 1) {
        nmi = 0.0;
    } else {
        ga.resize(qa);
        for(int q = 0; q < qa; q++) ga[q]=0;
        gb.resize(qb);
        for(int q = 0; q < qb; q++) gb[q]=0;

        vector< vector<int> > A;
        vector< vector<int> > B;
        A.resize(qa); // Existing structure
        B.resize(qa); // Counting structure
        for(int i = 0; i < n; i++) {
            int q = pa[i];
            int t = pb[i];
            ga[q]++;
            gb[t]++;
            int idx = -1;
            for(int j = 0; j < A[q].size(); j++) {
                if(A[q][j] == t) {
                    idx=j;
                    break;
                }
            }
            if(idx == -1) { // Pair [x y] did not show up
                A[q].push_back(t);
                B[q].push_back(1);
            } else { // [x y] is there
                B[q][idx] += 1;
            }
        }
        double Ha = 0;
        for(int q = 0; q < qa; q++) {
            if(ga[q] == 0) continue;
            double prob = 1.0*ga[q]/n;
            Ha += prob*log(prob);
        }
        double Hb=0;
        for(int q = 0; q < qb; q++) {
            if(gb[q] == 0) continue;
            double prob = 1.0*gb[q]/n;
            Hb += prob*log(prob);
        }
        double Iab=0;
        for(int q = 0; q < qa; q++) {
            for(int idx = 0; idx < A[q].size(); idx++) {
                double prob = 1.0*B[q][idx]/n;
                int t = A[q][idx];
                Iab += prob*log(prob/ ( 1.0*ga[q]/n*gb[t]/n ));
            }
        }

        nmi = -2.0*Iab/(Ha+Hb);
    }
}

void Solution::ComputeCentroidIndex(Solution* ground_truth) {
	int d = pb_data.GetD();
	int m = pb_data.GetM();
	double dist, mindist;
	int cmin = -1;
	vector<bool> orphan(m, true);
	double ci = m;

    for(int i = 0; i < m; i++) {
		mindist = MAX_FLOAT;
		for(int j = 0; j < m; j++) {
			dist = SquaredEuclidean(centroids[i], ground_truth->GetCentroids(j), d);
			if(dist < mindist) {
				mindist = dist;
				cmin = j;
			}
		}
        if(orphan[cmin] == true) {
            ci--;
        }
		orphan[cmin] = false;
	}

	centroid_index = ci;
}

void Solution::ComputeExternalMetrics(unsigned short* assignment) {
    Solution* ground_truth = new Solution(assignment, 0.0, pb_data);
    ComputeCRand(ground_truth);
    ComputeNmi(ground_truth);
    ComputeCentroidIndex(ground_truth);
    delete ground_truth;
}