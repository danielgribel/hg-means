#include "GeneticOperations.h"

GeneticOperations::GeneticOperations(PbData pb_data, Param param) {
    this->pb_data = pb_data;
    this->param = param;
}

GeneticOperations::~GeneticOperations() {
    for(int i = 0; i < population.size(); i++) {
        DeleteSolution(i);
    }
}

Solution* GeneticOperations::SelectParent() {
    Solution* parent_solution = NULL;
    double best_cost = MAX_FLOAT;
    int r;
    for(int i = 0; i < param.w; i++) {
        r = rand() % population.size();
        if(GetCost(r) < best_cost) {
            parent_solution = population[r];
            best_cost = GetCost(r);
        }
    }
    return parent_solution;
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
    double best_cost = MAX_FLOAT;
    for(int i = 0; i < param.size_population; i++) {
        unsigned short* assignment = GetKmeansAssignment(x);    
        double alpha = 1.0;
        if(param.mutation) {
            alpha = RandBetween(0.0, 1.0);
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

vector<int> GeneticOperations::GetCardinality(int** cluster_size) {
    vector<int> cardinality (pb_data.GetM());
    for(int i = 0; i < pb_data.GetM(); i++) {
        cardinality[i] = cluster_size[0][i];
    }
    return cardinality;
}

void GeneticOperations::PushMax(vector< pair<double, int> >& heap, double cost, int val) {
    heap.push_back( pair<double, int>(cost, val) );
    push_heap(heap.begin(), heap.end());
}

int GeneticOperations::PopMax(vector< pair<double, int> >& heap) {
    make_heap(heap.begin(), heap.end());
    pop_heap(heap.begin(), heap.end());
    heap.pop_back();
    return heap.front().second;
}

pair<double, int> GeneticOperations::FrontMax(vector< pair<double, int> >& heap) {
    make_heap(heap.begin(), heap.end());
    return heap.front();
}

// O(max_population x n)
void GeneticOperations::SelectSurvivors(const Dataset* x) {
    int m = pb_data.GetM();
    unordered_set<Item, KeyHash, KeyEqual> hash_table;
    int max_population = population.size();
    vector<Solution*> new_population;
    vector<int> discarded(max_population);
    vector< pair<double, int> > heap_individuals;
    vector< pair<double, int> > heap_clones;
    Kmeans* algorithm = new HamerlyKmeans();

    for(int i = 0; i < max_population; i++) {
        algorithm->initialize(x, m, GetAssignment(i), 1); // O(m)
        vector<int> card = GetCardinality(algorithm->getClusterSize()); // O(m)
        sort(card.begin(), card.end());
        Item an_item;
        an_item.cost = GetCost(i);
        an_item.cardinality = card;
        if (hash_table.find(an_item) != hash_table.end()) {
            PushMax(heap_clones, GetCost(i), i);
        } else {
            hash_table.insert(an_item);
            PushMax(heap_individuals, GetCost(i), i);
        }
        discarded[i] = 0;
    }
    delete algorithm;
    
    int id;
    int j = 0;
    
    // For the two whiles: O(max_population - size_population)
    while((j < (max_population - param.size_population)) && (heap_clones.size() > 0)) {
        id = FrontMax(heap_clones).second;
        PopMax(heap_clones);
        discarded[id] = 1;
        j++;
    }

    while(j < (max_population - param.size_population)) {
        id = FrontMax(heap_individuals).second;
        PopMax(heap_individuals);
        discarded[id] = 1;
        j++;
    }

    for(int i = 0; i < max_population; i++) { // O(max_population)
        if(discarded[i] == 0) {
            new_population.push_back(population[i]);
        } else {
            delete population[i];
        }
    }
    population = new_population;
}

// Get the centroids assignment (perfect matching)
vector<long> GeneticOperations::MinAssignment(vector< vector<double> > c1, vector< vector<double> > c2) { //O(m^2d)
    int m = pb_data.GetM();
    int d = pb_data.GetD();
    dlib::matrix<double> cost(m,m);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            cost(i,j) = -1.0 * SquaredEuclidean(c1[i], c2[j], d);
        }
    }
    vector<long> assignment = max_cost_assignment(cost);
    return assignment;
}

Solution* GeneticOperations::Crossover(Solution* p1, Solution* p2) {
    int d = pb_data.GetD();
    int m = pb_data.GetM();
    double alpha = 0.5 * (p1->GetAlpha() + p2->GetAlpha());

    vector< vector<double> > c1 = p1->GetCentroids();
    vector< vector<double> > c2 = p2->GetCentroids();
    vector< vector<double> > c3 (m, vector<double>(d));

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
    Solution* offspring = new Solution(c3, alpha, pb_data);
    offspring->Repair();
    return offspring;
}

void GeneticOperations::HGMeans(const Dataset* x) {
    int it = 0;
    int last_improvement = 0;

    // Generate the initial population
    CreateInitialPopulation(x);

    // Genetic algorithm general loop
	while(((it-last_improvement) < param.no_improvement_it) && (it < param.max_it)) {

        // Select parent solutions for mate
        Solution* p1 = SelectParent();
        Solution* p2 = SelectParent();

        // Apply the crossover
        Solution* current_solution = Crossover(p1, p2);

        // Mutate mutation factor
        if(param.mutation) {
            current_solution->MutateAlpha();
        }

        // Apply the mutation mutation
        current_solution->Mutate();

        // Perform local search (K-means)
        current_solution->DoLocalSearch(x);

        // Add child solution to population
        AddSolution(current_solution);

        // Update the best solution
        if(current_solution->GetCost() < best_solution->GetCost()) {
            ReplaceBestSolution(current_solution);
            last_improvement = it;
        }

        // Select the survivors
        if(population.size() > param.max_population) {
            SelectSurvivors(x);
        }
        it++;
    }
}