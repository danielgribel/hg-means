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
        
        int DoHash(int* size, int m);
        
    public:
        
        Hash();
        
        void Insert(Item newItem, int m);
        
        bool Exist(int* cardinality, double cost, int m);

        bool CheckSize(int* cardinalityA, int* cardinalityB, int m);
        
        ~Hash();
};

#endif