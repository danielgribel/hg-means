/* */
#include "GeneticOperations.h"
#include "PbRun.h"
#include "Evaluator.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <map>
#include <ctime>
#include <sys/resource.h>
#include <unistd.h>
#include <algorithm>
#include <queue>

using namespace std;

#define SAVE_FILE true
#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

PbRun * gaLoop(Dataset const *x, unsigned short m, Param prm) {
    clock_t begin = clock();

    const int n = x->n;
    int it = 0;
    int lastImprovement = 0;
    double elapsedSecs;
    double nbIter = 0.0;
    double bestCost = MathUtils::MAX_FLOAT;
    unsigned short* bestSolution = new unsigned short[n];
    double alpha, bestAlpha;

    GeneticOperations* genetic = new GeneticOperations();

    // Generate the initial population
    vector<Solution*> population = genetic->GetInitialPopulation(x, m, prm);

    // Store the best solution
    for(unsigned short i = 0; i < population.size(); i++) {
        if(population[i]->GetCost() < bestCost) {
            bestCost = population[i]->GetCost();
            bestAlpha = population[i]->GetAlpha();
            copy(population[i]->GetAssignment(), population[i]->GetAssignment() + n, bestSolution);
        }
    }

    // Genetic algorithm general loop
	while(((it-lastImprovement) < prm.itNoImprovement) && (it < prm.maxIt)) {

        // Select parent solutions for mate
        Solution* p1 = genetic->SelectParent(population, prm.W);
        Solution* p2 = genetic->SelectParent(population, prm.W);

        // Apply the crossover
        Solution* current_solution = genetic->Crossover(p1, p2, x, m, 0.5 * (p1->GetAlpha() + p2->GetAlpha()));
        
        // Mutate mutation factor
        if(prm.mutation) {
            current_solution->MutateAlpha();
        }

        // Apply the mutation
        current_solution->Mutate();

        // Perform local search (K-means)
        current_solution->DoLocalSearch(x);

        // Add child solution to population
        population.push_back(current_solution);

        // Update the best solution if is the case
        if(current_solution->GetCost() < bestCost) {
            bestCost = current_solution->GetCost();
            bestAlpha = current_solution->GetAlpha();
            copy(current_solution->GetAssignment(), current_solution->GetAssignment() + n, bestSolution);
            lastImprovement = it;
        }

        // Select the survivors
        if(population.size() > prm.maxPopulation) {
            population = genetic->SelectSurvivors(population, prm.sizePopulation, x, m);
        }
        it++;
    }
    elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
    Solution * best = new Solution(bestSolution, bestCost, bestAlpha, x, m);
    PbRun * sol = new PbRun(best, elapsedSecs);
    
    for(unsigned short i = 0; i < population.size(); i++) {
        delete population[i];
    }

    cout << m << " " << bestCost << " " << elapsedSecs << endl;

    delete genetic;

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
        x->data[i] = 1.0*x->data[i];
    }

    if (x == NULL) {
        cerr << "Please load a dataset first" << endl;
        return;
    }

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
        PbRun * r = gaLoop(x, m, prm);

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
                Solution* y_pred = r->GetSolution();
                Solution* y_ = new Solution(y, 0.0, x, m);
                Evaluator * eval = new Evaluator(n, m, d,
                    y_pred->GetAssignment(), y, y_pred->GetCentroids(), y_->GetCentroids());                
                crand = eval->cRand();
                nmi = eval->nmi();
                ci = eval->centroidIndex();
                delete eval;
                delete y_;
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