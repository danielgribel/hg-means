#ifndef PbData_H
#define PbData_H

#include <iostream>

struct Param {
	int W;
	int sizePopulation;
    int maxPopulation;
    int maxIt;
    int itNoImprovement;
    int nbRuns;
    bool mutation;
    bool eval;
};

class PbData {
    
    private:

        double* data;

        int n;

        int d;

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