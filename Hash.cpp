#include "Hash.h"

Hash::Hash() {
    length = 37;
    array = new std::vector<Item>[length];
}

int Hash::hash(int* cardinality, int m) {
    std::sort(cardinality, cardinality + m);
    int value = 0;

    for(int i = 0; i < m; i++) {
        value = value + (i+1)*cardinality[i];
    }

    return value % length;
}

void Hash::insert(Item newItem, int m) {
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

bool Hash::exist(int* cardinality, double cost, int m) {
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

Hash::~Hash() {
    for(int i = 0; i < length; i++) {
        for(int j = 0; j < array[i].size(); j++) {
            delete [] array[i][j].cardinality;
        }
    }
    delete [] array;
}