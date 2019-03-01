/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 */

#include "HGMeans.h"

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
            genetic->GetBestSolution()->ComputeExternalMetrics();
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
    param.w = 2;
    param.mutation = 1;
    param.nb_runs = 1;
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

HGMeans::HGMeans() {

}

HGMeans::~HGMeans() {

}

void HGMeans::Go(char* filename, int size_population, int it, const std::vector<int>& m) {
    string filenme_str(filename);
    InputValidator command(filenme_str, size_population, it, m);

    if(command.Validate()) {
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
                Run(16007, pb_data, x);
            } else {
                cerr << "The number of clusters must be less than " << x->n << " (number of data points)" << endl;
            }
        }
        pb_data.DeleteGroundTruth();
        delete x;
    }
}

int main(int argc, char** argv) {
    std::vector<int> nb_clusters;
    if(argc >= 5) {
        for(int i = 4; i < argc; i++) {
            nb_clusters.push_back(atoi(argv[i]));
        }
        HGMeans hgmeans;
        hgmeans.Go(argv[1], atoi(argv[2]), atoi(argv[3]), nb_clusters);
    } else {
        cerr << "Insufficient number of parameters provided. Please use the following input format:" << 
        endl << "./hgmeans DatasetPath Pi_min N2 [M]" << endl;
    }
    return 0;
}