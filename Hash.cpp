#include "Hash.h"

Hash::Hash() {
    length = 37;
    array = vector< vector<Item> >(length);
}

int Hash::DoHash(int* cardinality, int m) {
    sort(cardinality, cardinality + m);
    int value = 0;

    for(int i = 0; i < m; i++) {
        value = value + (i+1)*cardinality[i];
    }

    return value % length;
}

void Hash::Insert(Item newItem, int m) {
    int idx = DoHash(newItem.cardinality, m);
    array[idx].push_back(newItem);
}

bool Hash::CheckSize(int* cardinalityA, int* cardinalityB, int m) {
    for(int i = 0; i < m; i++) {
        if(cardinalityA[i] != cardinalityB[i]) {
            return false;
        }
    }
    return true;
}

bool Hash::Exist(int* cardinality, double cost, int m) {
    int idx = DoHash(cardinality, m);
    for(int i = 0; i < array[idx].size(); i++) {
        if ((CheckSize(array[idx][i].cardinality, cardinality, m)) &&
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
}