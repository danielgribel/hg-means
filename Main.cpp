/* */
#include "GeneticOperations.h"
#include "Evaluator.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;

#define SAVE_FILE true
#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

void Run(int seed, string fileData, Param prm, int m) {
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

    ofstream writer_output;
    stringstream filename_output;

    fileData = ReplaceString(fileData, DATA_PATH, "");
    fileData = ReplaceString(fileData, ".txt", "");

    filename_output << "out/" << fileData << '_' <<
                setw(3) << setfill('0') << m << "_" <<
                prm.size_population << '_' <<
                prm.max_it << ".out";

    // Clean file content if it already exists
    if(SAVE_FILE) {
        writer_output.open (filename_output.str().c_str());
        writer_output.close();
    }

    for(int i = 0; i < prm.nb_runs; i++) {
        clock_t begin = clock();

        // Create GeneticOperations instance
        GeneticOperations* genetic = new GeneticOperations(pb_data, prm);

        // Run HG-Means algorithm
        genetic->HGMeans(x);
        
        // Measure cpu time in seconds
        double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
        
        cout << fileData << " " << pb_data.GetM() << " " << genetic->GetBestSolution()->GetCost() << " " << elapsedSecs << endl;
        
        if(SAVE_FILE) {
            writer_output.open (filename_output.str().c_str(), ofstream::out | ofstream::app);  
            writer_output << prm.size_population << " ";
            writer_output << prm.max_it << " ";
            writer_output << fileData << " ";
            writer_output << m << " ";
            writer_output << fixed << setprecision(10) << genetic->GetBestSolution()->GetCost() << " ";
            writer_output << fixed << setprecision(4) << elapsedSecs << " ";
        }
        if(prm.eval) {
            string fileLabels = CLASS_PATH + fileData + ".txt";
            // Open the labels file
            ifstream inputLabels(fileLabels.c_str());
            if (! inputLabels) {
                cerr << "Unable to open labels file: " << fileLabels << endl;
                if(SAVE_FILE) {
                    writer_output << "\n";
                    writer_output.close();
                }
                return;
            }
            int max = 0;
            unsigned short* y = new unsigned short[n];
            for (int i = 0; i < n; ++i) {
                inputLabels >> y[i];
                y[i] = y[i] - 1; // Labels file starts from 1
                if(y[i] > max)
                    max = y[i];
            }
            if(max != m-1) {
                cerr << "Number of labels does not match m" << endl;
                if(SAVE_FILE) {
                    writer_output << "\n";
                    writer_output.close();
                }
                return;
            } else {
                Solution* ground_truth = new Solution(y, 0.0, pb_data);
                Evaluator * eval = new Evaluator(pb_data, genetic->GetBestSolution(), ground_truth);                
                // Clustering indexes
                if(SAVE_FILE) {
                    writer_output << fixed << setprecision(4) << eval->CRand() << " ";
                    writer_output << fixed << setprecision(4) << eval->Nmi() << " ";
                    writer_output << fixed << setprecision(4) << eval->CentroidIndex();
                }
                delete eval;
                delete ground_truth;
            }
        }
        if(SAVE_FILE) {
            writer_output << "\n";
            writer_output.close();
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
        parameters.size_population = atoi(argv[2]);
        parameters.max_it = atoi(argv[3]);
        parameters.eval = atoi(argv[4]);
        parameters.w = 2;
        parameters.mutation = 1;
        parameters.nb_runs = 1;
        parameters.max_population = 2*parameters.size_population;
        parameters.no_improvement_it = parameters.max_it/10;
        if(parameters.no_improvement_it < 1) {
            parameters.no_improvement_it = 1;
        }
        for(int i = 5; i < argc; i++) {
            m = atoi(argv[i]);
            if(m > 0 && m < USHRT_MAX) { // Range for m should be [1, 65535]
                Run(16007, fileName, parameters, m);
            } else {
                cerr << "The number of clusters is out of the limit [1, " << USHRT_MAX << "]" << endl;
            }
        }
    } catch (const exception& e) {
        cerr << e.what();
    }
    return 0;
}