#ifndef PbData_H
#define PbData_H

#include <iostream>
#include <string>

using namespace std;

// Parameters of the problem
struct Param {

    // Size of Tournament selection for parents mate
	int w;

    // Population size
	int size_population;

    // Maximum size allowed to population
    int max_population;

    // Maximum number of iteratios the algorithm will take
    int max_it;

    // Maximum number of iteratios without improvement the algorithm will take
    int no_improvement_it;

    // Number of runs
    int nb_runs;

    // If mutation is activated or not
    bool mutation;

    // If evaluation is activated or not -- in terms of C-Rand, NMI, Centroid Index, etc
    bool eval;
};

class PbData {
    
    private:

        // Instance (dataset) name
        string instance_name;

        // Features vector: Represented as a linearized vector of size n x d
        double* data;

        // Number of data points (samples)
        int n;

        // Number of dimensions (features)
        int d;

        // Number of clusters
        int m;
    
    public:

        PbData(string instance_name, double* data, int n, int d, int m) {
            this->instance_name = instance_name;
            this->data = data;
            this->n = n;
            this->d = d;
            this->m = m;
        };

        PbData() {};

        ~PbData() {};

        string GetInstanceName() { return instance_name; };

        double* GetData() { return data; };

        int GetN() { return n; };

        int GetD() { return d; };

        int GetM() { return m; };
};

#endif