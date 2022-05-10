/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 */

#include "HGMeans.h"

void SaveHeader(ofstream& writer_output, stringstream& filename_output, PbData pb_data) {
    writer_output.open (filename_output.str().c_str());
    // Write the header
    string header = string("DATASET ") +
                    "NB_CLUSTERS " +
                    "SIZE_POPULATION " +
                    "MAX_IT " +
                    "OBJECTIVE " +
                    "TIME";

    if(pb_data.GetParam().eval && pb_data.GetNbClasses() == pb_data.GetM()) {
        header = header + " CRAND " + "NMI " + "CI";
    }
    
    for (int i = 0; i < pb_data.GetN(); ++i) {
        header = header + " X" + to_string(i+1);
    }
    
    writer_output << header << " " << "\n";
    writer_output.close();
}

void SaveOutput(ofstream& writer_output, stringstream& filename_output, GeneticOperations* genetic, double elapsedSecs) {
    Param param = genetic->GetParam();
    int m = genetic->GetPbData().GetM();
    string instance = genetic->GetPbData().GetInstanceName();
    Solution* solution = genetic->GetBestSolution();
    
    writer_output.open (filename_output.str().c_str(), ofstream::out | ofstream::app);  

    writer_output << instance << " ";
    writer_output << m << " ";
    writer_output << param.size_population << " ";
    writer_output << param.max_it << " ";
    writer_output << fixed << setprecision(10) << solution->GetCost() << " ";
    writer_output << fixed << setprecision(4) << elapsedSecs;

    // Attach clustering indexes to output
    if(param.eval && genetic->GetPbData().GetNbClasses() == m) {
        writer_output << " ";
        writer_output << fixed << setprecision(4) << solution->GetCRand() << " ";
        writer_output << fixed << setprecision(4) << solution->GetNmi() << " ";
        writer_output << fixed << setprecision(4) << solution->GetCentroidIndex();
    }

    for (int i = 0; i < genetic->GetPbData().GetN(); ++i) {
        writer_output << " " << solution->GetAssignment()[i] + 1;
    }

    writer_output << "\n";
    writer_output.close();
}

void PrintResult(GeneticOperations* genetic, double cpu_time) {
    Solution* best_solution = genetic->GetBestSolution();
    PbData pb_data = genetic->GetPbData();
    Param param = genetic->GetParam();
    
    std::cout << std::fixed;
    cout << "-- Optimization finished." << endl;
    cout << "   Solution objective: " << setprecision(4) << best_solution->GetCost() << endl
         << "   CPU time (s): " << setprecision(2) << cpu_time << endl;
    
    if(param.eval && pb_data.GetNbClasses() == pb_data.GetM()) {
        cout << "   Clustering performance: "
             << setprecision(4) << best_solution->GetCRand() << " (C-Rand), "
             << setprecision(4) << best_solution->GetNmi() << " (NMI), "
             << setprecision(4) << best_solution->GetCentroidIndex() << " (CI)" << endl;
    }
    cout << endl;
}

string FilenameOutput(PbData pb_data, bool save) {
    stringstream filename_output;

    string out_folder = "hgm_out/";
    mode_t mode = 0733; // UNIX style permissions
    int result = 0;

    if(save) {
        #if defined(_WIN32)
            result = _mkdir(out_folder.c_str()); // can be used on Windows
        #else 
            result = mkdir(out_folder.c_str(), mode); // can be used on non-Windows
        #endif
        // if (result != 0) {
        //     cout << "No output directory created at /hg-means. It already exists, or writing permissions should be checked." << endl;
        // }
    }

    filename_output << out_folder << pb_data.GetInstanceName() << '_' <<
                setw(3) << setfill('0') << pb_data.GetM() << "_" <<
                pb_data.GetParam().size_population << '_' <<
                pb_data.GetParam().max_it << ".out";
    
    return filename_output.str();
}

void RunFile(int seed, PbData pb_data, const Dataset* x, bool save) {
    srand(seed);

    ofstream writer_output;
    stringstream filename_output;
    filename_output << FilenameOutput(pb_data, save);
    
    // Clean file content if it already exists
    if(save) {
        SaveHeader(writer_output, filename_output, pb_data);
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
            genetic->GetBestSolution()->ComputeExternalMetrics();
        }

        PrintResult(genetic, cpu_time);

        if(save) {
            // Save output results 
            SaveOutput(writer_output, filename_output, genetic, cpu_time);
        }
        delete genetic;
    }
}

std::vector<int> Run(int seed, PbData pb_data, const Dataset* x, bool save) {
    srand(seed);

    ofstream writer_output;
    stringstream filename_output;
    filename_output << FilenameOutput(pb_data, save);
    
    // Clean file content if it already exists
    if(save) {
        SaveHeader(writer_output, filename_output, pb_data);
    }
    
    std::vector<int> outcome(x->n);
    double best_obj = MAX_FLOAT;

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
            genetic->GetBestSolution()->ComputeExternalMetrics();
        }

        PrintResult(genetic, cpu_time);

        if(save) {
            // Save output results 
            SaveOutput(writer_output, filename_output, genetic, cpu_time);
        }

        Solution* solution = genetic->GetBestSolution();
        double obj = solution->GetCost();
        if(obj < best_obj) {
            best_obj = obj;
            for(int i = 0; i < x->n; i++) {
                outcome[i] = solution->GetAssignment()[i];
            }
        }

        delete genetic;
    }

    return outcome;
}

Param LoadParam(InputValidator command) {
    Param param;
    param.size_population = command.GetPiMin();
    param.max_it = command.GetMaxIt();
    param.w = 2;
    param.mutation = 1;
    param.nb_runs = command.GetNbIt();
    param.eval = true;
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

const Dataset* LoadData(std::vector< std::vector<double> > dataset) {
    int n = dataset.size();
    int d = dataset[0].size();

    // Allocate storage
    const Dataset *x = new Dataset(n, d);
    
    // Read the data values directly into the dataset
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            x->data[d*(i % n) + j] = dataset[i][j];
        }
    }

    return x;
}

HGMeans::HGMeans() {

}

HGMeans::~HGMeans() {

}

void HGMeans::GoFile(char* filename, int size_population, int max_it, int nb_it, const std::vector<int>& m, bool save) {
    string filenme_str(filename);
    InputValidator command(filenme_str, size_population, max_it, nb_it, m);

    if(command.Validate(true)) {
        Param param = LoadParam(command);
        const Dataset* x = LoadData(command.GetDatasetPath());
        if (x == NULL) {
            cerr << "Please load a dataset first" << endl;
            return;
        }
        PbData pb_data(command.GetLabelPath(), command.GetDatasetName(), x->data, x->n, x->d, param);
        for(int i = 0; i < command.GetNbClusters().size(); i++) {
            if(command.GetNbClusters()[i] < pb_data.GetN()) { // m <= n
                pb_data.SetM(command.GetNbClusters()[i]);
                RunFile(16007, pb_data, x, save);
            } else {
                cerr << "The number of clusters must be less than " << x->n << " (number of data points)" << endl;
            }
        }
        pb_data.DeleteGroundTruth();
        delete x;
    }
}

std::vector< std::vector<int> > HGMeans::Go(std::vector< std::vector<double> > dataset, std::vector<int> y, int size_population, int max_it, int nb_it, const std::vector<int>& m, bool save) {
    string filenme_str("temp");
    InputValidator command(filenme_str, size_population, max_it, nb_it, m);

    std::vector< std::vector<int> > results;

    if(command.Validate(false)) {
        Param param = LoadParam(command);
        const Dataset* x = LoadData(dataset);
        if (x == NULL) {
            cerr << "Please load a dataset first" << endl;
            return { {NULL} };
        }
        PbData pb_data(y, command.GetDatasetName(), x->data, x->n, x->d, param);
        for(int i = 0; i < command.GetNbClusters().size(); i++) {
            if(command.GetNbClusters()[i] < pb_data.GetN()) { // m <= n
                pb_data.SetM(command.GetNbClusters()[i]);
                std::vector<int> res = Run(16007, pb_data, x, save);
                results.push_back(res);
            } else {
                cerr << "The number of clusters must be less than " << x->n << " (number of data points)" << endl;
            }
        }

        // std::cout << "The X:" << "\n";
        // for (int i = 0; i < x->n; ++i) {
        //     for (int j = 0; j < x->d; ++j) {
        //         std::cout << x->data[x->d*(i % x->n) + j] << " ";
        //     }
        //     std::cout << "\n";
        // }

        // std::cout << "The Y:" << "\n";
        // for (int i = 0; i < x->n; ++i) {
        //     std::cout << pb_data.GetTruthAssignment()[i] << " ";
        // }
        // std::cout << "\n";

        pb_data.DeleteGroundTruth();
        delete x;
    }

    return results;
}

int main(int argc, char** argv) {
    std::vector<int> nb_clusters;

    string w = argv[argc - 1];

    bool save = false;
    int lim_clusters = argc;

    if(w == "w") {
        save = true;
        lim_clusters -= 1;
    }

    if(argc >= 6) {
        for(int i = 5; i < lim_clusters; i++) {
            nb_clusters.push_back(atoi(argv[i]));
        }
        HGMeans hgmeans;
        hgmeans.GoFile(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), nb_clusters, save);
    } else {
        cerr << "Insufficient number of parameters provided. Please use the following input format:" << 
        endl << "./hgmeans DATASET_PATH PI_MIN N2 NB_IT [M]" << endl;
    }
    return 0;
}