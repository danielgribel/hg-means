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
                if (dataset_path.find(DATA_PATH) != std::string::npos &&
                    dataset_path.find(".txt") != std::string::npos) {
                    dataset_name = ReplaceString(dataset_path, DATA_PATH, "");
                    dataset_name = ReplaceString(dataset_name, ".txt", "");
                }
            }
        }

        string GetDatasetPath() { return dataset_path; };

        string GetDatasetName() { return dataset_name; };

        int GetPiMin() { return pi_min; };

        int GetMaxIt() { return max_it; };

        bool IsExternalEval() { return external_eval; };

        vector<int> GetNbClusters() { return nb_clusters; };

        bool Validate() {

            if(nb_arg < 6) {
                cerr << "Insufficient number of parameters" << endl;
                return false;
            }
            
            ifstream input(dataset_path.c_str());

            if (!input) {
                cerr << "Unable to open data file: " << dataset_path << endl;
                return false;
            }

            // Defining bounds on variables
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
    if(prm.eval && genetic->GetPbData().GetNbClasses() == m) {
        writer_output << " ";
        writer_output << fixed << setprecision(4) << solution->GetCRand() << " ";
        writer_output << fixed << setprecision(4) << solution->GetNmi() << " ";
        writer_output << fixed << setprecision(4) << solution->GetCentroidIndex();
    }
    writer_output << "\n";
    writer_output.close();
}

void Run(int seed, PbData pb_data, const Dataset* x) {
    srand(seed);

    ofstream writer_output;
    stringstream filename_output;

    filename_output << "out/" << pb_data.GetInstanceName() << '_' <<
                setw(3) << setfill('0') << pb_data.GetM() << "_" <<
                pb_data.GetParam().size_population << '_' <<
                pb_data.GetParam().max_it << ".out";

    // Clean file content if it already exists
    if(SAVE_FILE) {
        writer_output.open (filename_output.str().c_str());
        writer_output.close();
    }

    for(int i = 0; i < pb_data.GetParam().nb_runs; i++) {
        clock_t begin = clock();

        // Create GeneticOperations instance
        GeneticOperations* genetic = new GeneticOperations(pb_data);

        // Run HG-Means algorithm
        genetic->HGMeans(x);
        
        Solution * s = genetic->GetBestSolution();

        // Measure cpu time in seconds
        double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;

        cout << pb_data.GetInstanceName() << " " << pb_data.GetM() << " " << s->GetCost() << " " << elapsedSecs;
        
        if (pb_data.GetParam().eval && pb_data.GetM() == pb_data.GetNbClasses()) {
            s->ComputeExternalMetrics(pb_data.GetTruthAssignment());
            cout << " " << s->GetCRand() << " " << s->GetNmi() << " " << s->GetCentroidIndex();
        }
        cout << endl;

        if(SAVE_FILE) {
            // Save output results 
            SaveOutput(writer_output, filename_output, genetic, elapsedSecs);
        }
        delete genetic;
    }
}

Param LoadParam(InputValidator command) {
    Param param;
    param.size_population = command.GetPiMin();
    param.max_it = command.GetMaxIt();
    param.eval = command.IsExternalEval();
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
            pb_data.SetM(command.GetNbClusters()[i]);
            Run(16007, pb_data, x);
        }
        delete x;
    }
    return 0;
}