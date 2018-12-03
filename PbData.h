#ifndef PbData_H
#define PbData_H

#include <iostream>

// Parameters of the problem
struct Param {
	int w;
	int size_population;
    int max_population;
    int max_it;
    int no_improvement_it;
    int nb_runs;
    bool mutation;
    bool eval;
};

class PbData {
    
    private:

        // Features vector: Represented as a linearized vector of size n x d
        double* data;

        // Number of data points (samples)
        int n;

        // Number of dimensions (features)
        int d;

        // Number of clusters
        int m;
    
    public:

        PbData(double* data, int n, int d, int m) {
            this->data = data;
            this->n = n;
            this->d = d;
            this->m = m;
        };

        PbData() {};

        ~PbData() {};

        double* GetData() { return data; };

        int GetN() { return n; };

        int GetD() { return d; };

        int GetM() { return m; };
};

#endif