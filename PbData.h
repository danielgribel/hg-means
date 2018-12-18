#ifndef PbData_H
#define PbData_H

/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 */

#include <iostream>
#include <string>
#include <fstream>

#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

using namespace std;

// Stores the parameters of the problem
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

/* PbData stores the information of the problem: instance name, features vector,
 * number of data points, number of dimensions, and number of clusters
 */
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

        Param param;

        unsigned short* truth_assignment;

        int nb_classes;
    
    public:

        PbData(string instance_name, double* data, int n, int d, Param param) {
            this->instance_name = instance_name;
            this->data = data;
            this->n = n;
            this->d = d;
            this->param = param;
            if(param.eval) {
                LoadGroundTruth();
            }
        };

        PbData() {};

        ~PbData() {};

        string GetInstanceName() { return instance_name; };

        double* GetData() { return data; };

        int GetN() { return n; };

        int GetD() { return d; };

        int GetM() { return m; };

        int GetNbClasses() { return nb_classes; };

        Param GetParam() { return param; };

        unsigned short* GetTruthAssignment() { return truth_assignment; };

        void SetM(int m) { this->m = m; };

        void LoadGroundTruth() {
            string file_labels = CLASS_PATH + instance_name + ".txt";
            ifstream inputLabels(file_labels.c_str());
            if (! inputLabels) { // Verify if labels file was successfully loaded
                cerr << "Unable to open labels file: " << file_labels << endl;
                param.eval = false;
            } else {
                nb_classes = 0;
                truth_assignment = new unsigned short[n];
                for (int i = 0; i < n; ++i) {
                    inputLabels >> truth_assignment[i];
                    truth_assignment[i] = truth_assignment[i] - 1; // Labels file starts from 1
                    if(truth_assignment[i] > nb_classes)
                        nb_classes = truth_assignment[i];
                }
                nb_classes++;
            }
        }
};

#endif