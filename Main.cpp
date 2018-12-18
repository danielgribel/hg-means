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

class InputValidator {
    
    private:

        int nb_arg;

        string dataset_path;

        string dataset_name;

        int pi_min;

        int max_it;

        bool external_eval;

        vector<int> nb_clusters;

    public:

        InputValidator(int argc, char* argv[]) {

            nb_arg = argc; 

            if(nb_arg >= 6) {
                // Dataset path
                dataset_path = argv[1];

                // Population size (Pi_min paramter)
                pi_min = atoi(argv[2]);

                // Maxixum number of iterations parameter
                max_it = atoi(argv[3]);

                // External clustering evaluation parameter
                external_eval = atoi(argv[4]);

                for(int i = 5; i < argc; i++) {
                    nb_clusters.push_back(atoi(argv[i]));
                }
                if (dataset_path.find(DATA_PATH) != string::npos &&
                    dataset_path.find(".txt") != string::npos) {
                    dataset_name = ReplaceString(dataset_path, DATA_PATH, "");
                    dataset_name = ReplaceString(dataset_name, ".txt", "");
                }
            }
        }

        string GetDatasetPath() { return dataset_path; };

        string GetDatasetName() { return dataset_name; };

        int GetPiMin() { return pi_min; };

        int GetMaxIt() { return max_it; };

        bool ExternalEval() { return external_eval; };

        vector<int> GetNbClusters() { return nb_clusters; };

        bool Validate() {

            // At least 6 parameters should be provided:
            // Program name, dataset, size of population, maximum iterations, evaluation and number of clusters
            if(nb_arg < 6) {
                cerr << "Insufficient number of parameters provided. Please use the following input format:" << 
                endl << "./hgmeans DatasetPath Pi_min N2 Evaluate [M]" << endl;
                return false;
            }
            
            ifstream input(dataset_path.c_str());

            if (!input) {
                cerr << "Unable to open data file: " << dataset_path << endl;
                return false;
            }

            // Checking bounds of variables
            if(pi_min < 1 || pi_min > 100000) {
                cerr << "Pi_min out of bounds" << endl;
                return false;
            }
            
            if(max_it < 1 || max_it > 100000) {
                cerr << "Max_it out of bounds" << endl;
                return false;
            }

            if(external_eval != 0 && external_eval != 1) {
                cerr << "Incorrect value for evaluation" << endl;
                return false;
            }
                
            for(int i = 0; i < nb_clusters.size(); i++) {
                if(nb_clusters[i] < 1 || nb_clusters[i] > USHRT_MAX) { // Range for m must be [1, 65535]
                    cerr << "The number of clusters is out of the limit [1, " << USHRT_MAX << "]" << endl;
                    return false;
                }
            }
            return true;
        }
};

void SaveOutput(ofstream& writer_output, stringstream& filename_output, GeneticOperations* genetic, double elapsedSecs) {
    Param param = genetic->GetParam();
    int m = genetic->GetPbData().GetM();
    string instance = genetic->GetPbData().GetInstanceName();
    Solution* solution = genetic->GetBestSolution();
    
    writer_output.open (filename_output.str().c_str(), ofstream::out | ofstream::app);  
    writer_output << param.size_population << " ";
    writer_output << param.max_it << " ";
    writer_output << instance << " ";
    writer_output << m << " ";
    writer_output << fixed << setprecision(10) << solution->GetCost() << " ";
    writer_output << fixed << setprecision(4) << elapsedSecs;

    // Attach clustering indexes to output
    if(param.eval && genetic->GetPbData().GetNbClasses() == m) {
        writer_output << " ";
        writer_output << fixed << setprecision(4) << solution->GetCRand() << " ";
        writer_output << fixed << setprecision(4) << solution->GetNmi() << " ";
        writer_output << fixed << setprecision(4) << solution->GetCentroidIndex();
    }
    writer_output << "\n";
    writer_output.close();
}

void PrintResult(GeneticOperations* genetic, double cpu_time) {
    Solution* best_solution = genetic->GetBestSolution();
    PbData pb_data = genetic->GetPbData();
    Param param = genetic->GetParam();

    cout << "-- Optimization finished: Solution objective = " << best_solution->GetCost() << " | CPU time = " << cpu_time << " s";
    
    if(param.eval && pb_data.GetNbClasses() == pb_data.GetM()) {
        cout << setprecision(4) << " | C-Rand = " << best_solution->GetCRand();
        cout << setprecision(4) << " | NMI = " << best_solution->GetNmi();
        cout << setprecision(4) << " | CI = " << best_solution->GetCentroidIndex();
    }
    cout << endl << endl;
}

string FilenameOutput(PbData pb_data) {
    stringstream filename_output;

    filename_output << "out/" << pb_data.GetInstanceName() << '_' <<
                setw(3) << setfill('0') << pb_data.GetM() << "_" <<
                pb_data.GetParam().size_population << '_' <<
                pb_data.GetParam().max_it << ".out";
    
    return filename_output.str();
}

void Run(int seed, PbData pb_data, const Dataset* x) {
    srand(seed);

    ofstream writer_output;
    stringstream filename_output;
    filename_output << FilenameOutput(pb_data);
    
    // Clean file content if it already exists
    if(SAVE_FILE) {
        writer_output.open (filename_output.str().c_str());
        writer_output.close();
    }

    for(int i = 0; i < pb_data.GetParam().nb_runs; i++) {
        clock_t begin = clock();

        // Create GeneticOperations instance
        GeneticOperations* genetic = new GeneticOperations(pb_data);

        cout << "-- Starting optimization: " << pb_data.GetInstanceName()  << " dataset | m = " << pb_data.GetM() << " clusters" << endl;

        // Run HG-Means algorithm
        genetic->HGMeans(x);
        
        // Measure cpu time in seconds
        double cpu_time = double(clock() - begin) / CLOCKS_PER_SEC;

        if (pb_data.GetParam().eval && pb_data.GetM() == pb_data.GetNbClasses()) {
            genetic->GetBestSolution()->ComputeExternalMetrics(pb_data.GetTruthAssignment());
        }

        PrintResult(genetic, cpu_time);

        if(SAVE_FILE) {
            // Save output results 
            SaveOutput(writer_output, filename_output, genetic, cpu_time);
        }
        delete genetic;
    }
}

Param LoadParam(InputValidator command) {
    Param param;
    param.size_population = command.GetPiMin();
    param.max_it = command.GetMaxIt();
    param.eval = command.ExternalEval();
    param.w = 2;
    param.mutation = 1;
    param.nb_runs = 1;
    param.max_population = 2*param.size_population;
    param.no_improvement_it = param.max_it/10;
    if(param.no_improvement_it < 1) {
        param.no_improvement_it = 1;
    }
    return param;
}

const Dataset* LoadData(string data_path) {
    ifstream input(data_path.c_str());

    // Read number of points and dimensionality of data
    int n, d;
    input >> n;
    input >> d;

    // Allocate storage
    const Dataset *x = new Dataset(n, d);
    
    // Read the data values directly into the dataset
    for (int i = 0; i < n * d; ++i) {
        input >> x->data[i];
    }

    return x;
}

int main(int argc, char** argv) {

    InputValidator command(argc, argv);

    if(command.Validate()) {
        
        Param param = LoadParam(command);

        const Dataset* x = LoadData(command.GetDatasetPath());

        if (x == NULL) {
            cerr << "Please load a dataset first" << endl;
            return 0;
        }

        PbData pb_data(command.GetDatasetName(), x->data, x->n, x->d, param);

        for(int i = 0; i < command.GetNbClusters().size(); i++) {
            if(command.GetNbClusters()[i] <= pb_data.GetN()) { // m <= n
                pb_data.SetM(command.GetNbClusters()[i]);
                Run(16007, pb_data, x);
            }
        }
        delete x;
    }
    return 0;
}