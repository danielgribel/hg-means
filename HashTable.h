//*****************************************************************
//  HashTable.h
//  HashTable
//
//  Created by Karlina Beringer on June 18, 2014.
//
//  This header file contains the Hash Table class declaration.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#ifndef HashTable_h
#define HashTable_h

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#define EPS 0.00000001

struct Item {
    double cost;
    int* cardinality;
};

//*****************************************************************
// Hash Table objects store a fixed number of Linked Lists.
//*****************************************************************
class HashTable {
    
    private:
        
        // Array is a reference to an array of Linked Lists.
        std::vector<Item>* array;
        
        // Length is the size of the Hash Table array.
        int tableLength;
        
        // Returns an array location for a given item key.
        int hash(int* size, int m);
        
    public:
        
        // Constructs the empty Hash Table object.
        // Array length is set to 37 by default.
        HashTable();
        
        // Adds an item to the Hash Table.
        void insertItem(Item newItem, int m);
        
        // Returns true if the operation is successful.
        bool existItem(int* cardinality, double cost, int m);
        
        void destroy();
        
        // De-allocates all memory used for the Hash Table.
        ~HashTable();
};

#endif