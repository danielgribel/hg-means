/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 */

#include "GeneticOperations.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;

#define SAVE_FILE true
#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

void SaveOutput(ofstream& writer_output, stringstream& filename_output, GeneticOperations* genetic, Solution* ground_truth, double elapsedSecs) {
    Param prm = genetic->GetParam();
    int m = genetic->GetPbData().GetM();
    string instance = genetic->GetPbData().GetInstanceName();
    Solution* solution = genetic->GetBestSolution();
    
    writer_output.open (filename_output.str().c_str(), ofstream::out | ofstream::app);  
    writer_output << prm.size_population << " ";
    writer_output << prm.max_it << " ";
    writer_output << instance << " ";
    writer_output << m << " ";
    writer_output << fixed << setprecision(10) << solution->GetCost() << " ";
    writer_output << fixed << setprecision(4) << elapsedSecs;

    // Attach clustering indexes to output
    if(prm.eval) {
        writer_output << " ";
        writer_output << fixed << setprecision(4) << solution->GetCRand() << " ";
        writer_output << fixed << setprecision(4) << solution->GetNmi() << " ";
        writer_output << fixed << setprecision(4) << solution->GetCentroidIndex();
    }
    writer_output << "\n";
    writer_output.close();
}

void PrintResult(Solution* best_solution, Solution* ground_truth, double cpu_time, bool eval) {
    cout << "-- Optimization finished: Solution objective = " << best_solution->GetCost() << " | CPU time = " << cpu_time << " s";
    if(eval) {
        best_solution->ComputeExternalMetrics(ground_truth);
        cout << setprecision(4) << " | C-Rand = " << best_solution->GetCRand();
        cout << setprecision(4) << " | NMI = " << best_solution->GetNmi();
        cout << setprecision(4) << " | CI = " << best_solution->GetCentroidIndex();
    }
    cout << endl << endl;
}

void Run(int seed, string file_data, Param prm, int m) {
    srand(seed);

    // Open the data file
    ifstream input(file_data.c_str());
    if (! input) {
        cerr << "Unable to open data file: " << file_data << endl;
        return;
    }
    
    // Read number of points and dimensionality of data
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

    file_data = ReplaceString(file_data, DATA_PATH, "");
    file_data = ReplaceString(file_data, ".txt", "");

    // Create an instance of PbData
    PbData pb_data(file_data, x->data, x->n, x->d, m);

    ofstream writer_output;
    stringstream filename_output;

    filename_output << "out/" << file_data << '_' <<
                setw(3) << setfill('0') << m << "_" <<
                prm.size_population << '_' <<
                prm.max_it << ".out";

    // Clean file content if it already exists
    if(SAVE_FILE) {
        writer_output.open (filename_output.str().c_str());
        writer_output.close();
    }

    // Try to obtain ground truth labels -- open file, load labels, etc
    // If successful, create a ground truth solution; else, set prm.eval=False
    string instance = pb_data.GetInstanceName();
    Solution* ground_truth;

    if(prm.eval) {
        string file_labels = CLASS_PATH + instance + ".txt";
        ifstream inputLabels(file_labels.c_str());
        if (! inputLabels) { // Verify if labels file was successfully loaded
            cerr << "Unable to open labels file: " << file_labels << endl;
            prm.eval = false;
        } else {
            int max = 0;
            unsigned short* y = new unsigned short[n];
            for (int i = 0; i < n; ++i) {
                inputLabels >> y[i];
                y[i] = y[i] - 1; // Labels file starts from 1
                if(y[i] > max)
                    max = y[i];
            }
            if(max != m-1) { // Verify if labels file is valid
                cerr << "Number of labels does not match m" << endl;
                prm.eval = false;
                delete [] y;
            } else {
                ground_truth = new Solution(y, 0.0, pb_data);
            }
        }
    }

    for(int i = 0; i < prm.nb_runs; i++) {
        clock_t begin = clock();

        // Create GeneticOperations instance
        GeneticOperations* genetic = new GeneticOperations(pb_data, prm);

        cout << "-- Starting optimization: " << file_data  << " dataset | m = " << pb_data.GetM() << " clusters" << endl;
        
        // Run HG-Means algorithm
        genetic->HGMeans(x);
        
        // Measure cpu time in seconds
        double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
        
        PrintResult(genetic->GetBestSolution(), ground_truth, elapsedSecs, prm.eval);
        
        if(SAVE_FILE) {
            // Save output results 
            SaveOutput(writer_output, filename_output, genetic, ground_truth, elapsedSecs);
        }
        delete genetic;
    }
    if(prm.eval) {
        delete ground_truth;
    }
    delete x;
}

int main(int argc, char** argv) {

    if(argc < 6) {
        cerr << "Insufficient number of parameters provided. Please use the following input format:" << 
        endl << "./hgmeans DatasetPath Pi_min N2 Evaluate [M]" << endl;
        return 0;
    }

    try {
        Param parameters;
        int m;
        string fileName = argv[1];
        parameters.size_population = atoi(argv[2]);
        parameters.max_it = atoi(argv[3]);
        parameters.eval = atoi(argv[4]);
        
        // Defining bounds on variables
        if(parameters.size_population < 1 || parameters.size_population > 100000) {
            parameters.size_population = 10; // Default
        }
        if(parameters.max_it < 1 || parameters.max_it > 100000) {
            parameters.max_it = 5000; // Default
        }
        if(parameters.eval != 0 && parameters.eval != 1) {
            parameters.eval = 0; // Default
        }

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
            if(m > 0 && m < USHRT_MAX) { // Range for m must be [1, 65535]
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