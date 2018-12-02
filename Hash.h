#ifndef Hash_h
#define Hash_h

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

#define EPS 0.00000001

struct Item {
    double cost;
    int* cardinality;
};

class Hash {
    
    private:

        vector< vector<Item> > array;
        
        int length;
        
        int hash(int* size, int m);
        
    public:
        
        Hash();
        
        void insert(Item newItem, int m);
        
        bool exist(int* cardinality, double cost, int m);
        
        ~Hash();
};

#endif