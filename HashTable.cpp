//*****************************************************************
//  HashTable.cpp
//  HashTable
//
//  Created by Kar Beringer on June 18, 2014.
//
//  This header file contains the Hash Table class definition.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#include "HashTable.h"

// Constructs the empty Hash Table object.
HashTable::HashTable() {
    tableLength = 37;
    array = new std::vector<Item>[tableLength];
}

// Returns an array location for a given item key.
int HashTable::hash(int* cardinality, int m) {
    std::sort(cardinality, cardinality + m);
    int value = 0;

    for(int i = 0; i < m; i++) {
        value = value + (i+1)*cardinality[i];
    }

    return value % tableLength;
}

// Adds an item to the Hash Table.
void HashTable::insertItem(Item newItem, int m) {
    int idx = hash(newItem.cardinality, m);
    array[idx].push_back(newItem);
}

bool checkSize(int* cardinalityA, int* cardinalityB, int m) {
    for(int i = 0; i < m; i++) {
        if(cardinalityA[i] != cardinalityB[i]) {
            return false;
        }
    }
    return true;
}

bool HashTable::existItem(int* cardinality, double cost, int m) {
    int idx = hash(cardinality, m);
    for(int i = 0; i < array[idx].size(); i++) {
        if ((checkSize(array[idx][i].cardinality, cardinality, m)) &&
            (cost > array[idx][i].cost - EPS) &&
            (cost < array[idx][i].cost + EPS)) {
            return true;
        }
    }
    return false;
}

void HashTable::destroy() {
    for(int i = 0; i < tableLength; i++) {
        for(int j = 0; j < array[i].size(); j++) {
            delete [] array[i][j].cardinality;
        }
    }
}

// De-allocates all memory used for the Hash Table.
HashTable::~HashTable() {
    destroy();
    delete [] array;
}