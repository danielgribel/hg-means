#include "GeneticOperations.h"

GeneticOperations::GeneticOperations(int size_population, int max_population, int w) {
    this->size_population = size_population;
    this->max_population = max_population;
    this->w = w;
}

GeneticOperations::~GeneticOperations() {
    for(unsigned short i = 0; i < population.size(); i++) {
        DeleteSolution(i);
    }
}

Solution* GeneticOperations::SelectParent() {
    Solution* bestSolution = NULL;
    double bestCost = MathUtils::MAX_FLOAT;
    unsigned short r;
    for(int i = 0; i < w; i++) {
        r = rand() % population.size();
        if(GetCost(r) < bestCost) {
            bestSolution = population[r];
            bestCost = GetCost(r);
        }
    }
    return bestSolution;
}

unsigned short* GeneticOperations::GetKmeansAssignment(const Dataset* x, unsigned short m) {
    Dataset *c = NULL;
    // K-means initialization
    c = init_centers(*x, m);
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;
    return assignment;
}

unsigned short* GeneticOperations::GetKppAssignment(const Dataset* x, unsigned short m) {
    Dataset *c = NULL;
    // K-means++ initialization
    c = init_centers_kmeanspp_v2(*x, m);
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;
    return assignment;
}

void GeneticOperations::CreateInitialPopulation(const Dataset* x, unsigned short m, bool is_mutable) {
    for(int i = 0; i < size_population; i++) {
        unsigned short* assignment = GetKmeansAssignment(x, m);    
        double alpha = 1.0;
        if(is_mutable) {
            alpha = MathUtils::fRand(0.0, 1.0);
        }
        Solution* s = new Solution(assignment, alpha, x, m);
        s->DoLocalSearch(x);
        AddSolution(s);
    }
}

int* GeneticOperations::GetCardinality(int** clusterSize, int m) {
    int* cardinality = new int[m];
    for(int i = 0; i < m; i++) {
        cardinality[i] = clusterSize[0][i];
    }
    return cardinality;
}

// O(Max_Pop x n)
void GeneticOperations::SelectSurvivors(Dataset const *x, int m) {
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
        algorithm->initialize(x, m, GetAssignment(i), 1); // O(m)
        int* card = GetCardinality(algorithm->getClusterSize(), m); // O(m)
        // Check if element is in hash: O(n) worst case
        if(table->exist(card, GetCost(i), m)) {
            heapClones->push_max(GetCost(i), i);
            delete [] card;
        } else {
            Item anItem;
            anItem.cost = GetCost(i);
            anItem.cardinality = card;
            table->insert(anItem, m);
            heapInd->push_max(GetCost(i), i);
        }
        discarded[i] = 0;
    }
    delete algorithm;
    
    int j = 0;
    
    // For the two whiles: O(maxPop - sizePop)
    while((j < (maxPopulation-size_population)) && (heapClones->getHeap().size() > 0)) {
        id = heapClones->front_max().second;
        heapClones->pop_max();
        discarded[id] = 1;
        j++;
    }

    while(j < (maxPopulation - size_population)) {
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

    population = newPopulation;
}

// Get the centroids assignment (perfect matching)
vector<long> GeneticOperations::MinAssignment(double** c1, double** c2, int m, int d) { //O(m^2d)
    dlib::matrix<double> cost(m,m);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            cost(i,j) = -1.0 * MathUtils::squaredEuclidean(c1[i], c2[j], d);
        }
    }
    std::vector<long> assignment = max_cost_assignment(cost);
    return assignment;
}

Solution* GeneticOperations::Crossover(Solution* p1, Solution* p2, const Dataset* x, const int m) {
    int d = x->d;    
    double alpha = 0.5 * (p1->GetAlpha() + p2->GetAlpha());
    double** c1 = p1->GetCentroids();
    double** c2 = p2->GetCentroids();
    double** c3 = new double* [m];
    
    for(int i = 0; i < m; i++) {
        c3[i] = new double[d];
    }

    vector<long> matching = MinAssignment(c1, c2, m, d); // O(m^3)  Hungarian method

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