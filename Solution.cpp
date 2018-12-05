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

Solution::Solution(vector< vector<double> > centroids, double cost, double alpha, PbData pb_data) {
	this->centroids = centroids;
    this->cost = cost;
    this->alpha = alpha;
    this->pb_data = pb_data;
    CentroidsToAssignment();
}

Solution::Solution(vector< vector<double> > centroids, double alpha, PbData pb_data) {
	this->centroids = centroids;
    this->alpha = alpha;
    this->pb_data = pb_data;
    CentroidsToAssignment();
}

Solution::~Solution() {
    delete [] assignment;
}

void Solution::AssignmentToCentroids() { // O(nd)
    InitCentroids();
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
    InitAssignment();
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();
    double dist, mindist;
    for(int i = 0; i < n; i++) {
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

    vector<bool> populated(m, false); // the (boolean) list of empty clusters 
    vector<int> empty_clusters; // the list of empty clusters
    
    // Identify populated and empty clusters
    for(int i = 0; i < n; i++)
        populated[assignment[i]] = true;

    // Construct the list of empty clusters
    for(int j = 0; j < m; j++) {
        if(populated[j] == false)
            empty_clusters.push_back(j);
    }

    if(empty_clusters.size() > 0) {
        vector<double> dist_centroid(n); // The distance from data points to their centroids
        double sum_dist = 0.0; // Sum of all distances from data points to their centroids
        vector<int> sizes(m, 0); // The size (cardinality) of clusters
        AssignmentToCentroids();
        vector<double> pr(n);

        for(int i = 0; i < n; i++) {
            dist_centroid[i] = PointCenterDist(i, centroids[assignment[i]], d, data);
            sum_dist = sum_dist + dist_centroid[i];
            sizes[assignment[i]] = sizes[assignment[i]]+1;
        }

        double z = 0.0;
        for(int i = 0; i < n; i++) {
            pr[i] = z + Pr(dist_centroid[i], sum_dist, alpha, n);
            z = pr[i];
        }

        int it = 0;
        int p;

        while(it < empty_clusters.size()) {
            double r = RandBetween(0.0, pr[n-1]); // O(1)
            int p = FindIndex(pr, r, 0, n-1) + 1;
            if(sizes[assignment[p]] > 1) {
                sizes[assignment[p]] = sizes[assignment[p]]-1;
                assignment[p] = empty_clusters[it];
                it++;
            }
        }
    }
	AssignmentToCentroids();
}

void Solution::Mutate() {
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double* data = pb_data.GetData();
    vector<double> dist_centroid (n);
    vector<double> pr (n);
    vector<int> barycenter_obj;
    vector< vector<double> > newcentroids (m, vector<double>(d));
    double dist, mindist;
    
    // Randomly select one centroid to remove it from the solution
    int barycenter = rand() % m;
    vector< vector<double> > c3 = centroids;

    // Keep the data points assigned to the centroid to be removed
    for(int i = 0; i < n; i++) { // O(n)
        if(assignment[i] == barycenter) {
            barycenter_obj.push_back(i);
        }
    }

    int b;
    // Re-assign data points to the closest remaining centroid
    for(int q = 0; q < barycenter_obj.size(); q++) { // O(|c|md)
        b = barycenter_obj[q];
        mindist = MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            if(j != barycenter) {
                dist = PointCenterDist(b, c3[j], d, data);
                if(dist < mindist) {
                    mindist = dist;
                    assignment[b] = j;
                }
            }
        }
    }

    double sumDist = 0.0;
    for(int i = 0; i < n; i++) {
        dist_centroid[i] = PointCenterDist(i, c3[assignment[i]], d, data);
        sumDist = sumDist + dist_centroid[i];
    }

    // Get the cumulative distance from each data point $i$ to its centroid defined in $mutation$, i.e., the solution with $m-1$ centroids
    double z = 0.0;
    for(int i = 0; i < n; i++) { // O(nd)
        pr[i] = z + Pr(dist_centroid[i], sumDist, alpha, n);
        z = pr[i];
    }

    // Wheel roulette random choose of a data point. data points far from their centroids are more likely to be chosen
    double r = RandBetween(0.0, pr[n-1]); // O(1)
    int p = FindIndex(pr, r, 0, n-1) + 1;

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
        if(PointCenterDist(i, newcentroids[barycenter], d, data) < dist_centroid[i]) {
            assignment[i] = barycenter;
        }
    }
    Repair();
}

void Solution::MutateAlpha() {
    alpha = alpha + RandBetween(-MUTATION_RATE, MUTATION_RATE);
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
    algorithm->initialize(x, pb_data.GetM(), assignment, 1);
    int numIt = algorithm->run(MAX_INT);
    cost = algorithm->getSSE();
    AssignmentToCentroids();
    delete algorithm;
}

void Solution::InitAssignment() {
	assignment = new unsigned short[pb_data.GetN()];
}

void Solution::InitCentroids() {
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    centroids = vector< vector<double> > (m, vector<double>(d, 0.0));
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

double Solution::Rand(Solution* ground_truth) {
	int a, b, c, d;
	CountRandCoefficients(ground_truth, a, b, c, d);
	double randIndex = 1.0*(a + d)/(a + b + c + d);
	return randIndex;
}

double Solution::CRand(Solution* ground_truth) {
	int a, b, c, d;
	CountRandCoefficients(ground_truth, a, b, c, d);
	int total = a + b + c + d;
	double crandIndex = (a - (1.0*(b + a)*(c + a))/total)/((1.0*(b + a + c + a))/2 - (1.0*(b + a)*(c + a))/total);
	return crandIndex;
}

// Implemented by @Carlo Nicolini
// More information in the original repository: https://github.com/CarloNicolini/rnmi
double Solution::Nmi(Solution* ground_truth) {
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
	if(qa == 1 && qb == 1) return 0.0;
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
	return -2.0*Iab/(Ha+Hb);
}

double Solution::CentroidIndex(Solution* ground_truth) {
	int d = pb_data.GetD();
	int m = pb_data.GetM();
	double dist, cmin;
	vector<int> orphan(m, 1);
	
	for(int i = 0; i < m; i++) {
		double mindist = MAX_FLOAT;
		for(int j = 0; j < m; j++) {
			dist = SquaredEuclidean(centroids[i], ground_truth->GetCentroids(j), d);
			if(dist < mindist) {
				mindist = dist;
				cmin = j;
			}
		}
		orphan[cmin] = 0;
	}

	double ci = 0.0;

	for(int i = 0; i < m; i++) {
		ci = ci + orphan[i];
	}

	return ci;
}