#include "GeneticOperations.h"

GeneticOperations::GeneticOperations(PbData pb_data, Param param) {
    this->pb_data = pb_data;
    this->param = param;
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
    for(int i = 0; i < param.W; i++) {
        r = rand() % population.size();
        if(GetCost(r) < bestCost) {
            bestSolution = population[r];
            bestCost = GetCost(r);
        }
    }
    return bestSolution;
}

unsigned short* GeneticOperations::GetKmeansAssignment(const Dataset* x) {
    Dataset *c = NULL;
    // K-means initialization
    c = init_centers(*x, pb_data.GetM());
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;
    return assignment;
}

unsigned short* GeneticOperations::GetKppAssignment(const Dataset* x) {
    Dataset *c = NULL;
    // K-means++ initialization
    c = init_centers_kmeanspp_v2(*x, pb_data.GetM());
    unsigned short* assignment = new unsigned short[x->n];
    assign(*x, *c, assignment);
    delete c;
    return assignment;
}

void GeneticOperations::CreateInitialPopulation(const Dataset* x) {
    Solution* best;
    double best_cost = MathUtils::MAX_FLOAT;
    for(int i = 0; i < param.sizePopulation; i++) {
        unsigned short* assignment = GetKmeansAssignment(x);    
        double alpha = 1.0;
        if(param.mutation) {
            alpha = MathUtils::RandBetween(0.0, 1.0);
        }
        Solution* s = new Solution(assignment, alpha, pb_data);
        s->DoLocalSearch(x);
        AddSolution(s);
        if(s->GetCost() < best_cost) {
            best = s;
            best_cost = s->GetCost();
        }
    }
    StoreBestSolution(best);
}

void GeneticOperations::StoreBestSolution(Solution* s) {
    unsigned short* best_assignment = new unsigned short[pb_data.GetN()];
    copy(s->GetAssignment(), s->GetAssignment() + pb_data.GetN(), best_assignment);
    best_solution = new Solution(best_assignment, s->GetCost(), s->GetAlpha(), pb_data);
}

void GeneticOperations::ReplaceBestSolution(Solution* s) {
    delete best_solution;
    StoreBestSolution(s);
}

int* GeneticOperations::GetCardinality(int** clusterSize) {
    int* cardinality = new int [pb_data.GetM()];
    for(int i = 0; i < pb_data.GetM(); i++) {
        cardinality[i] = clusterSize[0][i];
    }
    return cardinality;
}

// O(Max_Pop x n)
void GeneticOperations::SelectSurvivors(const Dataset* x) {
    int n = pb_data.GetN();
    int m = pb_data.GetM();
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
        int* card = GetCardinality(algorithm->getClusterSize()); // O(m)
        // Check if element is in hash: O(n) worst case
        if(table->Exist(card, GetCost(i), m)) {
            heapClones->PushMax(GetCost(i), i);
            delete [] card;
        } else {
            Item anItem;
            anItem.cost = GetCost(i);
            anItem.cardinality = card;
            table->Insert(anItem, m);
            heapInd->PushMax(GetCost(i), i);
        }
        discarded[i] = 0;
    }
    delete algorithm;
    
    int j = 0;
    
    // For the two whiles: O(maxPop - sizePop)
    while((j < (maxPopulation- param.sizePopulation)) && (heapClones->GetHeap().size() > 0)) {
        id = heapClones->FrontMax().second;
        heapClones->PopMax();
        discarded[id] = 1;
        j++;
    }

    while(j < (maxPopulation - param.sizePopulation)) {
        id = heapInd->FrontMax().second;
        heapInd->PopMax();
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
vector<long> GeneticOperations::MinAssignment(double** c1, double** c2) { //O(m^2d)
    int m = pb_data.GetM();
    int d = pb_data.GetD();
    dlib::matrix<double> cost(m,m);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            cost(i,j) = -1.0 * MathUtils::SquaredEuclidean(c1[i], c2[j], d);
        }
    }
    std::vector<long> assignment = max_cost_assignment(cost);
    return assignment;
}

Solution* GeneticOperations::Crossover(Solution* p1, Solution* p2) {
    int n = pb_data.GetN();
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double alpha = 0.5 * (p1->GetAlpha() + p2->GetAlpha());
    double** c1 = p1->GetCentroids();
    double** c2 = p2->GetCentroids();
    double** c3 = new double* [m];
    
    for(int i = 0; i < m; i++) {
        c3[i] = new double[d];
    }

    vector<long> matching = MinAssignment(c1, c2); // O(m^3)  Hungarian method

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

    Solution* off = new Solution(c3, alpha, pb_data);
    off->Repair();

    return off;
}

void GeneticOperations::HGMeans(const Dataset* x) {
    int it = 0;
    int lastImprovement = 0;

    // Generate the initial population
    CreateInitialPopulation(x);

    // Genetic algorithm general loop
	while(((it-lastImprovement) < param.itNoImprovement) && (it < param.maxIt)) {

        // Select parent solutions for mate
        Solution* p1 = SelectParent();
        Solution* p2 = SelectParent();

        // Apply the crossover
        Solution* current_solution = Crossover(p1, p2);
        
        // Mutate mutation factor
        if(param.mutation) {
            current_solution->MutateAlpha();
        }

        // Apply the mutation
        current_solution->Mutate();

        // Perform local search (K-means)
        current_solution->DoLocalSearch(x);

        // Add child solution to population
        AddSolution(current_solution);

        // Update the best solution
        if(current_solution->GetCost() < GetBestSolution()->GetCost()) {
            ReplaceBestSolution(current_solution);
            lastImprovement = it;
        }

        // Select the survivors
        if(GetPopulation().size() > param.maxPopulation) {
            SelectSurvivors(x);
        }
        it++;
    }
}