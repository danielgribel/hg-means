#include "Evaluator.h"
#include "MathUtils.h"

Evaluator::Evaluator(int n, int m, int d, unsigned short* y_pred, unsigned short* y, double** c_pred, double** c) {
	this->n = n;
	this->m = m;
	this->d = d;
	this->y_pred = y_pred;
	this->y = y;
	this->c_pred = c_pred;
	this->c = c;
}

Evaluator::~Evaluator() {

}

void Evaluator::countRandCoefficients(int& a, int& b, int& c, int& d) {
	a = 0;
	b = 0;
	c = 0;
	d = 0;

	for(int i = 0; i < n; i++) {
		for(int j = i+1; j < n; j++) {
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
double Evaluator::rand() {
	int a, b, c, d;
	countRandCoefficients(a, b, c, d);
	double randIndex = 1.0*(a + d)/(a + b + c + d);
	return randIndex;
}

// Get C-rand indicator (or adjusted rand), a measure for partitions agreement
double Evaluator::cRand() {
	int a, b, c, d;
	countRandCoefficients(a, b, c, d);
	int total = a + b + c + d;
	double crandIndex = (a - (1.0*(b + a)*(c + a))/total)/((1.0*(b + a + c + a))/2 - (1.0*(b + a)*(c + a))/total);
	return crandIndex;
}

// Get the Normalized mutual information indicator
double Evaluator::nmi() {
	unsigned short* pa = y_pred;
	unsigned short* pb = y;
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
	for(int i=0;i<n;i++){
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
double Evaluator::centroidIndex() {
	double dist, cmin;
	vector<int> orphan(m, 1);

	for(int i = 0; i < m; i++) {
		double mindist = MathUtils::MAX_FLOAT;
		for(int j = 0; j < m; j++) {
			dist = MathUtils::squaredEuclidean(c_pred[i], c[j], d);
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