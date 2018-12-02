#include "Evaluator.h"
#include "MathUtils.h"

Evaluator::Evaluator(PbData pb_data, Solution* solution, Solution* ground_truth) {
	this->pb_data = pb_data;
	this->solution = solution;
	this->ground_truth = ground_truth;
}

Evaluator::~Evaluator() {

}

void Evaluator::CountRandCoefficients(int& a, int& b, int& c, int& d) {
	a = 0;
	b = 0;
	c = 0;
	d = 0;

	unsigned short* y_pred = solution->GetAssignment();
	unsigned short* y = ground_truth->GetAssignment();

	for(int i = 0; i < pb_data.GetN(); i++) {
		for(int j = i+1; j < pb_data.GetN(); j++) {
			if( (y_pred[i] == y_pred[j]) && (y[i] == y[j]) ) {
				a++;
			}
			if( (y_pred[i] == y_pred[j]) && (y[i] != y[j]) ) {
				b++;
			}
			if( (y_pred[i] != y_pred[j]) && (y[i] == y[j]) ) {
				c++;
			}
			if( (y_pred[i] != y_pred[j]) && (y[i] != y[j]) ) {
				d++;
			}
		}
	}
}

// Get Rand indicator, a measure for partitions agreement
double Evaluator::Rand() {
	int a, b, c, d;
	CountRandCoefficients(a, b, c, d);
	double randIndex = 1.0*(a + d)/(a + b + c + d);
	return randIndex;
}

// Get C-rand indicator (or adjusted rand), a measure for partitions agreement
double Evaluator::CRand() {
	int a, b, c, d;
	CountRandCoefficients(a, b, c, d);
	int total = a + b + c + d;
	double crandIndex = (a - (1.0*(b + a)*(c + a))/total)/((1.0*(b + a + c + a))/2 - (1.0*(b + a)*(c + a))/total);
	return crandIndex;
}

// Get the Normalized mutual information indicator
double Evaluator::Nmi() {
	int n = pb_data.GetN();
	unsigned short* pa = solution->GetAssignment();
	unsigned short* pb = ground_truth->GetAssignment();
	int qa=-1,qb=-1;
	vector <int > ga;//group a
	vector <int > gb;//group b
	for(int i=0;i<n;i++){
		if(qa<pa[i]) qa=pa[i];
		if(qb<pb[i]) qb=pb[i];
	}
	qa++;
	qb++;
	if(qa==1 && qb==1) return 0.0;
	ga.resize(qa);
	for(int q=0;q<qa;q++) ga[q]=0;
	gb.resize(qb);
	for(int q=0;q<qb;q++) gb[q]=0;

	vector< vector<int> > A;
	vector< vector<int> > B;
	A.resize(qa); //existing structure
	B.resize(qa); //counting structure
	for(int i=0;i < n;i++){
		int q=pa[i];
		int t=pb[i];
		ga[q]++;
		gb[t]++;
		int idx=-1;
		for(int j=0;j<A[q].size();j++){
			if(A[q][j] == t) {
				idx=j;
				break;
			}
		}
		if(idx == -1){//pair [x y] did not show up
			A[q].push_back(t);
			B[q].push_back(1);
		}else{// [x y] is there
			B[q][idx] += 1;
		}
	}
	double Ha=0;
	for(int q=0;q<qa;q++){
		if(ga[q]==0) continue;
		double prob=1.0*ga[q]/n;
		Ha += prob*log(prob);
	}
	double Hb=0;
	for(int q=0;q<qb;q++){
		if(gb[q]==0) continue;
		double prob=1.0*gb[q]/n;
		Hb += prob*log(prob);
	}
	double Iab=0;
	for(int q=0;q<qa;q++){
		for(int idx=0;idx<A[q].size();idx++){
			double prob=1.0*B[q][idx]/n;
			int t=A[q][idx];
			Iab += prob*log(prob/ ( 1.0*ga[q]/n*gb[t]/n ));
		}
	}
	return -2.0*Iab/(Ha+Hb);
}

// Get the centroid index indicator
double Evaluator::CentroidIndex() {
	int d = pb_data.GetD();
	int m = pb_data.GetM();
	double dist, cmin;
	vector<int> orphan(m, 1);
	
	for(int i = 0; i < m; i++) {
		double mindist = MathUtils::MAX_FLOAT;
		for(int j = 0; j < m; j++) {
			dist = MathUtils::squaredEuclidean(solution->GetCentroids()[i], ground_truth->GetCentroids()[j], d);
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