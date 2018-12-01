#ifndef PbData_H
#define PbData_H

#include <iostream>

class PbData {
    
    private:

        double* data;

        int n;

        int d;

        int m;
    
    public:

        PbData();
        
        PbData(double* data, int n, int d, int m);

        ~PbData();

        double* GetData() { return data; };

        int GetN() { return n; };

        int GetD() { return d; };

        int GetM() { return m; };
};

#endif