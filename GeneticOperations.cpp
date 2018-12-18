#include "GeneticOperations.h"

GeneticOperations::GeneticOperations(PbData pb_data) {
    this->pb_data = pb_data;
    this->param = pb_data.GetParam();
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
        // Get random K-means initial solution
        unsigned short* assignment = GetKmeansAssignment(x);    
        
        double alpha = 1.0;

        // If mutation is enabled, choose alpha randomly
        if(param.mutation) {
            alpha = RandBetween(0.0, 1.0);
        }
        
        // Create new solution
        Solution* s = new Solution(assignment, alpha, pb_data);
        
        // Apply local search to solution
        s->DoLocalSearch(x);
        
        // Add solution to population
        AddSolution(s);

        if(s->GetCost() < best_cost) {
            best = s;
            best_cost = s->GetCost();
        }
    }
    // Save the best solution found so far
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

void GeneticOperations::DetectClones(const Dataset* x, vector< pair<double, int> >& heap_individuals, vector< pair<double, int> >& heap_clones) {
    
    // Hash table to detect clones
    unordered_set<Item, KeyHash, KeyEqual> hash_table;
    
    Kmeans* algorithm = new HamerlyKmeans();

    for(int i = 0; i < population.size(); i++) {
        
        algorithm->initialize(x, pb_data.GetM(), GetAssignment(i), 1); // O(m)
        
        // Get the cardinalities of clusters
        vector<int> card = GetCardinality(algorithm->getClusterSize()); // O(m)
        
        // Sort the list of cardinalities
        sort(card.begin(), card.end());
        
        // Create a hash item with the solution cost and the ordered cardinalities of clusters
        Item an_item;
        an_item.cost = GetCost(i);
        an_item.cardinality = card;

        if (hash_table.find(an_item) != hash_table.end()) {
            // If item is in the hash table, add it to the heap of clones
            PushMax(heap_clones, GetCost(i), i);
        } else {
            // If item is NOT in the hash table, add it to the hash table and to the heap of individuals
            hash_table.insert(an_item);
            PushMax(heap_individuals, GetCost(i), i);
        }
    }
    delete algorithm;
}

void GeneticOperations::ResetPopulation(vector< pair<double, int> >& heap_individuals, vector< pair<double, int> >& heap_clones) {
    
    // Current size of population
    int max_population = population.size();

    // The vector to be filled with the new population
    vector<Solution*> new_population;

    // Marks the solutions that will be discarded
    vector<int> discarded(max_population,0);

    // Counts the number of clone solutions
    int nb_clones = 0;

    // The two loops: O(max_population - size_population)
    // Remove clones until (max_population - size_population) solutions is achieved or no more clones exist
    for (int i = 0; (i < (max_population - param.size_population)) && (heap_clones.size() > 0); i++) {
        discarded[ FrontMax(heap_clones).second ] = 1;
        PopMax(heap_clones);
        nb_clones++;
    }

    // Remove non-clone solutions with bad-fitness
    for (int i = nb_clones; i < (max_population - param.size_population); i++) {
        discarded[ FrontMax(heap_individuals).second ] = 1;
        PopMax(heap_individuals);
    }

    // Build the new population
    for(int i = 0; i < max_population; i++) { // O(max_population)
        if(discarded[i] == 0) {
            new_population.push_back(population[i]);
        } else {
            delete population[i];
        }
    }
    population = new_population;
}

void GeneticOperations::SelectSurvivors(const Dataset* x) { // O(max_population x n)
    
    // Heap to detect individual solutions (non-clones)
    vector< pair<double, int> > heap_individuals;

    // Heap to detect clone solutions
    vector< pair<double, int> > heap_clones;
    
    // Detect clone and individual solutions and fill the heaps properly
    DetectClones(x, heap_individuals, heap_clones);

    // Remove clone/bad solutions and update the population
    ResetPopulation(heap_individuals, heap_clones);
}

vector<long> GeneticOperations::MinAssignment(vector< vector<double> > & c1, vector< vector<double> > & c2) {
    int m = pb_data.GetM();
    int d = pb_data.GetD();

    // Construct the matrix of distances between centroids
    dlib::matrix<double> cost(m,m);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            cost(i,j) = -1.0 * SquaredEuclidean(c1[i], c2[j], d);
        }
    }
    // Get the Min-Assignment from Dlib -- negative cost matrix is passed, so max_cost_assignment() is called
    vector<long> assignment = max_cost_assignment(cost);
    
    return assignment;
}

Solution* GeneticOperations::Crossover(Solution* p1, Solution* p2) {
    int m = pb_data.GetM();

    // Cross alpha parameter
    double alpha = 0.5 * (p1->GetAlpha() + p2->GetAlpha());

    // Get the centers of first parent
    vector< vector<double> > c1 = p1->GetCentroids();
    
    // Get the centers of second parent
    vector< vector<double> > c2 = p2->GetCentroids();
    
    // Initialize centers of child solution
    vector< vector<double> > c3 (m, vector<double>(pb_data.GetD()));

    // Hungarian method for Min-Assignment problem -- DLib implementation is called here
    vector<long> matching = MinAssignment(c1, c2); // O(m^3)

    // Create new centers -- for each assigned pair of centers, one is randomly chosen
    for(int i = 0; i < m; i++) { // O(md)
        if((rand() % 2) == 0) {
            c3[i] = c1[i];
        } else {
            c3[i] = c2[ matching[i] ];
        }
    }
    // Create child solution
    Solution* offspring = new Solution(c3, alpha, pb_data);

    // Repair solution in case of degenerated assignment 
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

        // Select the survivors if max population is achieved
        if(population.size() > param.max_population) {
            SelectSurvivors(x);
        }
        it++;
    }
}