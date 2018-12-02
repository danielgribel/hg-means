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

PbRun * ExecuteGA(Dataset const *x, PbData pb_data, Param prm) {
    clock_t begin = clock();

    const int n = pb_data.GetN();
    const int m = pb_data.GetM();
    int it = 0;
    int lastImprovement = 0;
    double elapsedSecs;

    // Create an instance of GeneticOperations
    GeneticOperations* genetic = new GeneticOperations(pb_data, prm.sizePopulation, prm.maxPopulation, prm.W);

    // Generate the initial population
    genetic->CreateInitialPopulation(x, prm.mutation);

    // Genetic algorithm general loop
	while(((it-lastImprovement) < prm.itNoImprovement) && (it < prm.maxIt)) {

        // Select parent solutions for mate
        Solution* p1 = genetic->SelectParent();
        Solution* p2 = genetic->SelectParent();

        // Apply the crossover
        Solution* current_solution = genetic->Crossover(p1, p2);
        
        // Mutate mutation factor
        if(prm.mutation) {
            current_solution->MutateAlpha();
        }

        // Apply the mutation
        current_solution->Mutate();

        // Perform local search (K-means)
        current_solution->DoLocalSearch(x);

        // Add child solution to population
        genetic->AddSolution(current_solution);

        // Update the best solution
        if(current_solution->GetCost() < genetic->GetBestSolution()->GetCost()) {
            genetic->ReplaceBestSolution(current_solution);
            lastImprovement = it;
        }

        // Select the survivors
        if(genetic->GetPopulation().size() > prm.maxPopulation) {
            genetic->SelectSurvivors(x);
        }
        it++;
    }
    elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
    PbRun * sol = new PbRun(genetic->GetBestSolution(), elapsedSecs);
    delete genetic;

    cout << m << " " << genetic->GetBestSolution()->GetCost() << " " << elapsedSecs << endl;

    return sol;
}

void demo(int seed, string fileData, Param prm, unsigned short m) {
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

    fileData = MathUtils::replaceString(fileData, DATA_PATH, "");
    fileData = MathUtils::replaceString(fileData, ".txt", "");

    outfile << "out/" << fileData << '_' <<
                setw(3) << setfill('0') << m << "_" <<
                prm.sizePopulation << '_' <<
                prm.maxIt << ".out";

    // Clean file content if it already exists
    if(SAVE_FILE) {
        myfile.open (outfile.str().c_str());
        myfile.close();
    }

    // Clustering indexes
    double crand = 0.0;
    double nmi = 0.0;
    double ci = 0.0;

    for(int i = 0; i < prm.nbRuns; i++) {
        PbRun * r = ExecuteGA(x, pb_data, prm);

        if(SAVE_FILE) {
            myfile.open (outfile.str().c_str(), ofstream::out | ofstream::app);  
            myfile << prm.sizePopulation << " ";
            myfile << prm.maxIt << " ";
            myfile << fileData << " ";
            myfile << m << " ";
            myfile << fixed << setprecision(10) << r->GetSolution()->GetCost() << " ";
            myfile << fixed << setprecision(4) << r->GetTime() << " ";
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
                Solution* solution = r->GetSolution();
                Solution* ground_truth = new Solution(y, 0.0, pb_data);
                Evaluator * eval = new Evaluator(pb_data, solution, ground_truth);                
                crand = eval->CRand();
                nmi = eval->Nmi();
                ci = eval->CentroidIndex();
                delete eval;
                delete ground_truth;
            }
        }

        if(SAVE_FILE) {
            if(prm.eval) {
                myfile << fixed << setprecision(4) << crand << " ";
                myfile << fixed << setprecision(4) << nmi << " ";
                myfile << fixed << setprecision(4) << ci;
            }
            myfile << "\n";
            myfile.close();
        }
        delete r->GetSolution();
        delete r;
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
            demo(16007, fileName, parameters, m);
        }

    } catch (const exception& e) {
        cerr << e.what();
    }
    return 0;
}