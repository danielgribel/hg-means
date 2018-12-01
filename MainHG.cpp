/* */

#include "hamerly/dataset.h"
#include "hamerly/general_functions.h"
#include "hamerly/hamerly_kmeans.h"
#include "Heap.h"
#include "Hash.h"
#include "Param.h"
#include "PbRun.h"
#include "Solution.h"
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
#include "dlib-master/dlib/optimization/max_cost_assignment.h"

using namespace std;

#define MUTATION_NOISE 0.2
#define SAVE_FILE true
#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

PbRun * gaLoop(Dataset const *x, unsigned short k, Param parameters);

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
            myfile << fixed << setprecision(10) << r->getSolution()->GetCost() << " ";
            myfile << fixed << setprecision(4) << r->getTime() << " ";
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
                Solution* y_pred = r->getSolution();
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
        delete r->getSolution();
        delete r;
    }

    delete x;
}

// Get the centroids assignment (perfect matching)
vector<long> minAssignment(double** c1, double** c2, int m, int d) { //O(m^2d)
    dlib::matrix<double> cost(m,m);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            cost(i,j) = -1.0 * MathUtils::squaredEuclidean(c1[i], c2[j], d);
        }
    }
    std::vector<long> assignment = max_cost_assignment(cost);
    return assignment;
}

Solution* crossover(Solution* p1, Solution* p2, const Dataset* x, const int m, double alpha) {
    int d = x->d;    
    double** c1 = p1->GetCentroids();
    double** c2 = p2->GetCentroids();
    double** c3 = new double* [m];
    
    for(int i = 0; i < m; i++) {
        c3[i] = new double[d];
    }

    vector<long> matching = minAssignment(c1, c2, m, d); // O(m^3)  Hungarian method

    for(int i = 0; i < m; i++) { // O(md)
        if((rand() % 2) == 0) {
            for(int j = 0; j < d; j++) {
                c3[i][j] = c1[i][j];
            }
        } else {
            for(int j = 0; j < d; j++) {
                c3[i][j] = c2[matching[i]][j];
            }
        }
    }

    Solution* off = new Solution(c3, alpha, x, m);
    off->Repair();

    return off;
}

unsigned short* getKmeansSolution(const Dataset* x, unsigned short m) {
    Dataset *c = NULL;

    // K-means initialization
    c = init_centers(*x, m);
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;

    return assignment;
}

unsigned short* getKppSolution(const Dataset* x, unsigned short m) {
    Dataset *c = NULL;

    // K-means++ initialization
    c = init_centers_kmeanspp_v2(*x, m);
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;
    
    return assignment;
}

vector<Solution*> getInitialPopulation(const Dataset* x, unsigned short m, Param prm) {
    vector<Solution*> population;
    for(int i = 0; i < prm.sizePopulation; i++) {
        unsigned short* assignment = getKmeansSolution(x, m);    
        double alpha = 1.0;
        if(prm.mutation) {
            alpha = MathUtils::fRand(0.0, 1.0);
        }
        Solution* s = new Solution(assignment, alpha, x, m);
        s->DoLocalSearch(x);
        population.push_back(s);
    }
    return population;
}

Solution* tournamentSelection(vector<Solution*> pop, int W) {
    Solution* bestSolution = NULL;
    double bestCost = MathUtils::MAX_FLOAT;
    unsigned short r;

    for(int i = 0; i < W; i++) {
        r = rand() % pop.size();
        if(pop[r]->GetCost() < bestCost) {
            bestSolution = pop[r];
            bestCost = pop[r]->GetCost();
        }
    }

    return bestSolution;
}

int* getCardinality(int** clusterSize, int m) {
    int* cardinality = new int[m];
    for(int i = 0; i < m; i++) {
        cardinality[i] = clusterSize[0][i];
    }
    return cardinality;
}

// O(Max_Pop x n)
vector<Solution*> selectSurvivors(vector<Solution*> population, const int sizePopulation, Dataset const *x, int m) {
    const int n = x->n;
    Hash* table = new Hash();
    int maxPopulation = population.size();
    vector<Solution*> newPopulation;
    vector<int> discarded(maxPopulation);
    int id;
    HeapPdi* heapInd = new HeapPdi();
    HeapPdi* heapClones = new HeapPdi(); 
    Kmeans* algorithm = new HamerlyKmeans();

    for(unsigned short i = 0; i < maxPopulation; i++) { // O(maxPop x n)
        algorithm->initialize(x, m, population[i]->GetAssignment(), 1); // O(m)
        int* card = getCardinality(algorithm->getClusterSize(), m); // O(m)
        // Check if element is in hash: O(n) worst case
        if(table->exist(card, population[i]->GetCost(), m)) {
            heapClones->push_max(population[i]->GetCost(), i);
            delete [] card;
        } else {
            Item anItem;
            anItem.cost = population[i]->GetCost();
            anItem.cardinality = card;
            table->insert(anItem, m);
            heapInd->push_max(population[i]->GetCost(), i);
        }
        discarded[i] = 0;
    }
    delete algorithm;
    
    int j = 0;
    
    // For the two whiles: O(maxPop - sizePop)
    while((j < (maxPopulation-sizePopulation)) && (heapClones->getHeap().size() > 0)) {
        id = heapClones->front_max().second;
        heapClones->pop_max();
        discarded[id] = 1;
        j++;
    }

    while(j < (maxPopulation-sizePopulation)) {
        id = heapInd->front_max().second;
        heapInd->pop_max();
        discarded[id] = 1;
        j++;
    }

    for(unsigned short i = 0; i < maxPopulation; i++) { // O(maxPop)
        if(discarded[i] == 0) {
            newPopulation.push_back(population[i]);
        } else {
            delete population[i];
        }
    }
    
    delete heapInd;
    delete heapClones;
    delete table;

    return newPopulation;
}

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

    // Generate the initial population
    vector<Solution*> population = getInitialPopulation(x, m, prm);
    
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
        Solution* p1 = tournamentSelection(population, prm.W);
        Solution* p2 = tournamentSelection(population, prm.W);

        // Apply the crossover
        Solution* current_solution = crossover(p1, p2, x, m, 0.5 * (p1->GetAlpha() + p2->GetAlpha()));
        
        // Mutate mutation factor
        if(prm.mutation) {
            current_solution->MutateAlpha();
        }

        // Apply the mutation
        current_solution->Mutate();

        // Perform local search (k-means)
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
            population = selectSurvivors(population, prm.sizePopulation, x, m);
        }

        it++;
    }

    elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
    nbIter = 1.0*nbIter/it;

    Solution * best = new Solution(bestSolution, bestCost, bestAlpha, x, m);
    PbRun * sol = new PbRun(best, elapsedSecs, lastImprovement, nbIter, 0.0, 1);
    
    for(unsigned short i = 0; i < population.size(); i++) {
        delete population[i];
    }

    cout << m << " " << bestCost << " " << elapsedSecs << endl;

    return sol;
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