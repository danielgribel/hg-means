/* */
#include "GeneticOperations.h"
#include "PbRun.h"
#include "Evaluator.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;

#define SAVE_FILE true
#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

void Run(int seed, string fileData, Param prm, unsigned short m) {
    srand(seed);

    // Open the data file
    ifstream input(fileData.c_str());
    if (! input) {
        cerr << "Unable to open data file: " << fileData << endl;
        return;
    }
    
    // Read the parameters
    int n, d;
    input >> n;
    input >> d;

    // Allocate storage
    Dataset *x = new Dataset(n, d);
    
    // Read the data values directly into the dataset
    for (int i = 0; i < n * d; ++i) {
        input >> x->data[i];
    }

    if (x == NULL) {
        cerr << "Please load a dataset first" << endl;
        return;
    }

    // Create an instance of PbData
    PbData pb_data(x->data, x->n, x->d, m);

    ofstream myfile;
    stringstream outfile;

    fileData = MathUtils::ReplaceString(fileData, DATA_PATH, "");
    fileData = MathUtils::ReplaceString(fileData, ".txt", "");

    outfile << "out/" << fileData << '_' <<
                setw(3) << setfill('0') << m << "_" <<
                prm.sizePopulation << '_' <<
                prm.maxIt << ".out";

    // Clean file content if it already exists
    if(SAVE_FILE) {
        myfile.open (outfile.str().c_str());
        myfile.close();
    }

    for(int i = 0; i < prm.nbRuns; i++) {
        clock_t begin = clock();

        // Create GeneticOperations instance
        GeneticOperations* genetic = new GeneticOperations(pb_data, prm);

        // Run the Genetic loop
        genetic->HGMeans(x);
        
        // Measure cpu time in seconds
        double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
        
        cout << pb_data.GetM() << " " << genetic->GetBestSolution()->GetCost() << " " << elapsedSecs << endl;
        
        if(SAVE_FILE) {
            myfile.open (outfile.str().c_str(), ofstream::out | ofstream::app);  
            myfile << prm.sizePopulation << " ";
            myfile << prm.maxIt << " ";
            myfile << fileData << " ";
            myfile << m << " ";
            myfile << fixed << setprecision(10) << genetic->GetBestSolution()->GetCost() << " ";
            myfile << fixed << setprecision(4) << elapsedSecs << " ";
        }

        if(prm.eval) {
            string fileLabels = CLASS_PATH + fileData + ".txt";
            // Open the labels file
            ifstream inputLabels(fileLabels.c_str());
            if (! inputLabels) {
                cerr << "Unable to open labels file: " << fileLabels << endl;
                if(SAVE_FILE) {
                    myfile << "\n";
                    myfile.close();
                }
                return;
            }

            int max = 0;
            unsigned short* y = new unsigned short[n];
            for (int i = 0; i < n; ++i) {
                inputLabels >> y[i];
                y[i] = y[i] - 1; // labels file starts from 1
                if(y[i] > max)
                    max = y[i];
            }

            if(max != m-1) {
                cerr << "Number of labels does not match m" << endl;
                if(SAVE_FILE) {
                    myfile << "\n";
                    myfile.close();
                }
                return;
            } else {
                Solution* ground_truth = new Solution(y, 0.0, pb_data);
                Evaluator * eval = new Evaluator(pb_data, genetic->GetBestSolution(), ground_truth);                
                // Clustering indexes
                if(SAVE_FILE) {
                    myfile << fixed << setprecision(4) << eval->CRand() << " ";
                    myfile << fixed << setprecision(4) << eval->Nmi() << " ";
                    myfile << fixed << setprecision(4) << eval->CentroidIndex();
                }
                delete eval;
                delete ground_truth;
            }
        }

        if(SAVE_FILE) {
            myfile << "\n";
            myfile.close();
        }
        delete genetic;
    }
    delete x;
}

int main(int argc, char** argv) {
    try {
        Param parameters;
        int m;
        string fileName = argv[1];
        parameters.sizePopulation = atoi(argv[2]);
        parameters.maxIt = atoi(argv[3]);
        parameters.eval = atoi(argv[4]);
        parameters.W = 2;
        parameters.mutation = 1;
        parameters.nbRuns = 1;
        parameters.maxPopulation = 2*parameters.sizePopulation;
        parameters.itNoImprovement = parameters.maxIt/10;
        if(parameters.itNoImprovement < 1) {
            parameters.itNoImprovement = 1;
        }
        for(int i = 5; i < argc; i++) {
            m = atoi(argv[i]);
            Run(16007, fileName, parameters, m);
        }

    } catch (const exception& e) {
        cerr << e.what();
    }
    return 0;
}