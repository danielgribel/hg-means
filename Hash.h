#ifndef Hash_h
#define Hash_h

/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 */

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

        // Hash table
        vector< vector<Item> > array;
        
        // Hash table length
        int length;
        
        // Apply hash function
        int DoHash(int* size, int m);
        
    public:
        
        Hash();
        
        // Insert a new item in the hash table
        void Insert(Item newItem, int m);
        
        // Check if the item exits in the table
        bool Exist(int* cardinality, double cost, int m);

        // Check if two solutions have the same clusters cardinality (cardinalities supposed to be sorted)
        bool CheckSize(int* cardinalityA, int* cardinalityB, int m);
        
        ~Hash();
};

#endif