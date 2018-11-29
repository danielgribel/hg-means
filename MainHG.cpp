/* Authors: Greg Hamerly and Jonathan Drake
 * Feedback: hamerly@cs.baylor.edu
 * See: http://cs.baylor.edu/~hamerly/software/kmeans.php
 * Copyright 2014
 */

/* This file is a common driver for running k-means implementations.
 *
 * The program's behavior is determined by a list of commands read from
 * standard input. Input lines should take one of the following forms:
 *
 * dataset some_dataset.txt
 * initialize k [random|kpp]
 * lloyd
 * hamerly
 * elkan
 * adaptive
 * annulus
 * compare
 * sort
 *
 * There are a number of shorthand alternatives,
 * e.g. init for initialize, data for dataset
 *
 * The dataset of n floating-point values in d-dimensional space is
 * read from the indicated file name and should have the form:
 *
 * n d
 * x(1,1) x(1,2) ... x(1, d)
 * .
 * .
 * .
 * x(n,1) x(n,2) ... x(n, d)
 *
 * Results go to standard output, and some extra
 * status information is also sent to standard error.
 */

#include "hamerly/dataset.h"
#include "hamerly/general_functions.h"
#include "hamerly/hamerly_kmeans.h"
#include "MinAssignment.h"
#include "Heap.h"
#include "HashTable.h"
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

using namespace std;

#define MUTATION_NOISE 0.2
#define SAVE_FILE true
#define DATA_PATH "data/"
#define CLASS_PATH "labels/"

PbRun * gaLoop(Dataset const *x, unsigned short k, Param parameters);
void fixSolution(unsigned short* solution, const Dataset* x, const int m, double alpha);

void deleteMatrix(double** matrix, int m) {
    for(int i = 0; i < m; i++) {
        delete [] matrix[i];
    }
    delete [] matrix;
}

int execute(Kmeans *algorithm, Dataset const *x, unsigned short k, unsigned short *assignment, int maxIterations) {

    int numThreads = 1;

    // Check for missing initialization
    if(assignment == NULL) {
        cerr << "Please initialize centers first" << endl;
        return 0;
    }
    if(x == NULL) {
        cerr << "Please load a dataset first" << endl;
        return 0;
    }

    // Time the execution and get the number of iterations
    algorithm->initialize(x, k, assignment, numThreads);
    int numIt = algorithm->run(maxIterations);

    return numIt;
}

double** assignmentToCentroids(unsigned short* solution, const Dataset* x, const int m) { // O(nd)
    const int n = x->n;
    const int d = x->d;
    double *data = x->data;
    double** centroid = new double*[m];
    vector<int> sizes(m);

    for(int i = 0; i < m; i++) {
        centroid[i] = new double[d];
        for(int j = 0; j < d; j++) {
            centroid[i][j] = 0.0;
        }
        sizes[i] = 0;
    }

    for(int i = 0; i < n; i++) {
        sizes[solution[i]] = sizes[solution[i]]+1;
        for(int j = 0; j < d; j++) {
            centroid[solution[i]][j] = centroid[solution[i]][j] + data[i*d+j];
        }
    }
    
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < d; j++) {
            if(sizes[i] != 0) {
                centroid[i][j] = (1.0*centroid[i][j])/sizes[i];
            }
        }
    }

    return centroid;
}

void centroidsToAssignment(unsigned short* solution, double** centroid, const Dataset* x, const int m) { // O(nmd)
    double mindist;
    double dist;

    for(int i = 0; i < x->n; i++) {
        mindist = MathUtils::MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            dist = MathUtils::pcDist(i, centroid[j], x->d, x->data);
            if(dist < mindist) {
                mindist = dist;
                solution[i] = j;
            }
        }
    }
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
            myfile << fixed << setprecision(10) << r->getSolution()->getCost() << " ";
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

                unsigned short* y_pred = r->getSolution()->getAssignment();
                double** c_pred = assignmentToCentroids(y_pred, x, m);
                double** c = assignmentToCentroids(y, x, m);

                Evaluator * eval = new Evaluator(n, m, d, y_pred, y, c_pred, c);

                crand = eval->cRand();
                nmi = eval->nmi();
                ci = eval->centroidIndex();

                deleteMatrix(c_pred, m);
                deleteMatrix(c, m);
                delete eval;
            }
            
            delete [] y;
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
double** assignment(double** c1, double** c2, int m, int d) { //O(m^2d)
    double** matrix = new double*[m];

    for(int i = 0; i < m; i++) {
        matrix[i] = new double[m];
    }

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            matrix[i][j] = MathUtils::squaredEuclidean(c1[i], c2[j], d);
        }
    }

    return matrix;
}

void crossover(unsigned short* offspring, unsigned short* p1, unsigned short* p2, const Dataset* x, const int m, double alpha) {
    int d = x->d;
    
    double** c1 = assignmentToCentroids(p1, x, m); // O(nd)
    double** c2 = assignmentToCentroids(p2, x, m); // O(nd)
    double** c3 = new double* [m];
    
    for(int i = 0; i < m; i++) {
        c3[i] = new double[d];
    }

    double** matrix = assignment(c1, c2, m, d); // O(m^2 d)
    vector<long> matching = minAssignment(matrix, m); // O(m^3)  Hungarian method

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

    centroidsToAssignment(offspring, c3, x, m); // O(nmd)
    fixSolution(offspring, x, m, alpha);

    deleteMatrix(c1, m);
    deleteMatrix(c2, m);
    deleteMatrix(c3, m);
    deleteMatrix(matrix, m);
}

void fixSolution(unsigned short* solution, const Dataset* x, const int m, double alpha) {
    vector<bool> populated(m, false); // the (boolean) list of empty clusters 
    vector<int> listEmpty; // the list of empty clusters
    
    // identify populated/empty clusters
    for(int i = 0; i < x->n; i++)
        populated[solution[i]] = true;

    // construct the list of empty clusters
    for(int j = 0; j < m; j++) {
        if(populated[j] == false)
            listEmpty.push_back(j);
    }

    if(listEmpty.size() > 0) {
        vector<double> distToCentroid(x->n); // the distance from data points to their centroids
        double sumDist = 0.0; // sum of all distances from data points to their centroids
        vector<int> sizes(m, 0); // the size (cardinality) of clusters
        double** centroids = assignmentToCentroids(solution, x, m);
        vector<double> pr(x->n);

        for(int i = 0; i < x->n; i++) {
            distToCentroid[i] = MathUtils::pcDist(i, centroids[solution[i]], x->d, x->data);
            sumDist = sumDist + distToCentroid[i];
            sizes[solution[i]] = sizes[solution[i]]+1;
        }

        double z = 0.0;
        for(int i = 0; i < x->n; i++) {
            pr[i] = z + MathUtils::Pr(distToCentroid[i], sumDist, alpha, x->n);
            z = pr[i];
        }

        int it = 0;
        int p;

        while(it < listEmpty.size()) {
            double r = MathUtils::fRand(0.0, pr[x->n-1]); // O(1)
            int p = MathUtils::findIndex(pr, r, 0, x->n-1) + 1;
            if(sizes[solution[p]] > 1) {
                sizes[solution[p]] = sizes[solution[p]]-1;
                solution[p] = listEmpty[it];
                it++;
            }
        }
        deleteMatrix(centroids, m);
    }
}

unsigned short* getKmeansSolution(const Dataset* x, unsigned short m) {
    Dataset *c = NULL;

    // k-means initialization
    c = init_centers(*x, m);
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;

    return assignment;
}

unsigned short* getKppSolution(const Dataset* x, unsigned short m) {
    Dataset *c = NULL;

    // k-means++ initialization
    c = init_centers_kmeanspp_v2(*x, m);
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;
    
    return assignment;
}

vector<Solution*> getInitialPopulation(const Dataset* x, unsigned short m, Param prm) {
    vector<Solution*> population;
    int numIt;
    Kmeans* algorithm = new HamerlyKmeans();

    for(int i = 0; i < prm.sizePopulation; i++) {
        unsigned short* assignment = getKmeansSolution(x, m);    
        numIt = execute(algorithm, x, m, assignment, MathUtils::MAX_INT);
        double alpha = 1.0;
        if(prm.mutation) {
            alpha = MathUtils::fRand(0.0, 1.0);
        }
        Solution* s = new Solution(assignment, algorithm->getSSE(), alpha);
        population.push_back(s);
    }
    delete algorithm;

    return population;
}

Solution* tournamentSelection(vector<Solution*> pop, int W) {
    Solution* bestSolution = NULL;
    double bestCost = MathUtils::MAX_FLOAT;
    unsigned short r;

    for(int i = 0; i < W; i++) {
        r = rand() % pop.size();
        if(pop[r]->getCost() < bestCost) {
            bestSolution = pop[r];
            bestCost = pop[r]->getCost();
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
    HashTable* table = new HashTable();
    int maxPopulation = population.size();
    vector<Solution*> newPopulation;
    vector<int> discarded(maxPopulation);
    int id;
    HeapPdi* heapInd = new HeapPdi();
    HeapPdi* heapClones = new HeapPdi(); 
    Kmeans* algorithm = new HamerlyKmeans();

    for(unsigned short i = 0; i < maxPopulation; i++) { // O(maxPop x n)
        algorithm->initialize(x, m, population[i]->getAssignment(), 1); // O(m)
        int* card = getCardinality(algorithm->getClusterSize(), m); // O(m)
        if(table->existItem(card, population[i]->getCost(), m)) { // check if element is in hash: O(n) worst case
            heapClones->push_max(population[i]->getCost(), i);
            delete [] card;
        } else {
            Item anItem;
            anItem.cost = population[i]->getCost();
            anItem.cardinality = card;
            table->insertItem(anItem, m);
            heapInd->push_max(population[i]->getCost(), i);
        }
        discarded[i] = 0;
    }
    delete algorithm;
    
    int j = 0;
    
    // for the two whiles: O(lamba), lambda = maxPop - sizePop
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

void mutation(unsigned short* offspring, unsigned short* mutated, const Dataset* x, const int m, double alpha) {
    vector<double> distCentroid(x->n);
    vector<double> pr (x->n);
    vector<int> barycenterObj;
    double** newcentroids = new double* [m];
    double mindist;
    double dist;
    int i;

    // randomly select one centroid to remove it from the solution
    int barycenter = rand() % m;

    double** c3 = assignmentToCentroids(offspring, x, m); // O(nd)

    // keep the data points assigned to the centroid to be removed
    for(int i = 0; i < x->n; i++) { // O(n)
        mutated[i] = offspring[i];
        if(offspring[i] == barycenter) {
            barycenterObj.push_back(i);
        }
    }

    // re-assign data points to the closest remaining centroid
    for(int q = 0; q < barycenterObj.size(); q++) { // O(|c|md)
        i = barycenterObj[q];
        mindist = MathUtils::MAX_FLOAT;
        for(int j = 0; j < m; j++) {
            if(j != barycenter) {
                dist = MathUtils::pcDist(i, c3[j], x->d, x->data);
                if(dist < mindist) {
                    mindist = dist;
                    mutated[i] = j;
                }    
            }
        }
    }

    double sumDist = 0.0;
    for(int i = 0; i < x->n; i++) {
        distCentroid[i] = MathUtils::pcDist(i, c3[mutated[i]], x->d, x->data);
        sumDist = sumDist + distCentroid[i];
    }

    // Get the cumulative distance from each data point $i$ to its centroid defined in $mutation$, i.e., the solution with $m-1$ centroids
    double z = 0.0;
    for(int i = 0; i < x->n; i++) { // O(nd)
        pr[i] = z + MathUtils::Pr(distCentroid[i], sumDist, alpha, x->n);
        z = pr[i];
    }

    // Wheel roulette random choose of a data point. data points far from their centroids are more likely to be chosen
    double r = MathUtils::fRand(0.0, pr[x->n-1]); // O(1)
    int p = MathUtils::findIndex(pr, r, 0, x->n-1) + 1;
    
    for(int i = 0; i < m; i++) {
        newcentroids[i] = new double[x->d];
    }

    // Re-insert the removed centroid in the position determined by the roulette wheel
    for(int i = 0; i < m; i++) { // O(md)
        if(i == barycenter) {
            for(int j = 0; j < x->d; j++) {
                newcentroids[barycenter][j] = x->data[(p*x->d)+j];
            }
        } else {
            for(int j = 0; j < x->d; j++) {
                newcentroids[i][j] = c3[i][j];
            }
        }
    }

    // Re-assign data points according to the solution with $m$ centroids
    for(int i = 0; i < x->n; i++) { // O(nd)
        if(MathUtils::pcDist(i, newcentroids[barycenter], x->d, x->data) < distCentroid[i]) {
            mutated[i] = barycenter;
        }
    }

    fixSolution(mutated, x, m, alpha);

    deleteMatrix(c3, m);
    deleteMatrix(newcentroids, m);
}

double mutationAlpha(double alpha) {
    alpha = alpha + MathUtils::fRand(-MUTATION_NOISE, MUTATION_NOISE);
    if(alpha > 1.0) {
        return 1.0;
    }
    if(alpha < 0.0) {
        return 0.0;
    }
    return alpha;
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
    double avgAlpha = 0.0;

    // Generate the initial population
    vector<Solution*> population = getInitialPopulation(x, m, prm);
    
    // Store the best solution
    for(unsigned short i = 0; i < population.size(); i++) {
        if(population[i]->getCost() < bestCost) {
            bestCost = population[i]->getCost();
            bestAlpha = population[i]->getAlpha();
            copy(population[i]->getAssignment(), population[i]->getAssignment() + n, bestSolution);
        }
    }

    Kmeans* algorithm = new HamerlyKmeans();

	while(((it-lastImprovement) < prm.itNoImprovement) && (it < prm.maxIt)) {

        // Select the parents for mate
        Solution* p1 = tournamentSelection(population, prm.W);
        Solution* p2 = tournamentSelection(population, prm.W);

        // Apply the crossover
        unsigned short* offspring = new unsigned short[n];
        alpha = 0.5 * (p1->getAlpha() + p2->getAlpha());
        crossover(offspring, p1->getAssignment(), p2->getAssignment(), x, m, alpha);
        
        // Apply the mutation
        unsigned short* mutated = new unsigned short[n];
        
        if(prm.mutation)
            alpha = mutationAlpha(alpha);

        mutation(offspring, mutated, x, m, alpha);
        
        avgAlpha = avgAlpha + alpha;

        // Apply local search to the mutated solution       
        nbIter = nbIter + execute(algorithm, x, m, mutated, MathUtils::MAX_INT);
        Solution* child = new Solution(mutated, algorithm->getSSE(), alpha);
        population.push_back(child);

        // Update the best solution if is the case
        if(child->getCost() < bestCost) {
            bestCost = child->getCost();
            bestAlpha = alpha;
            copy(child->getAssignment(), child->getAssignment() + n, bestSolution);
            lastImprovement = it;
        }

        // Select the survivors
        if(population.size() > prm.maxPopulation) {
            population = selectSurvivors(population, prm.sizePopulation, x, m);
        }

        it++;
        delete offspring;
    }
    delete algorithm;

    elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
    nbIter = 1.0*nbIter/it;
    avgAlpha = 1.0*avgAlpha/it;

    Solution * best = new Solution(bestSolution, bestCost, bestAlpha);
    PbRun * sol = new PbRun(best, elapsedSecs, lastImprovement, nbIter, avgAlpha, 1);
    
    for(unsigned short i = 0; i < population.size(); i++) {
        delete population[i];
    }

    cout << "CodeDaniel-Refac: " << m << " " << bestCost << " " << elapsedSecs << endl;

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