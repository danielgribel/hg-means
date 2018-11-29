//*****************************************************************
//  LinkedList.cpp
//  HashTable
//
//  Created by Karlina Beringer on June 16, 2014.
//
//  This header file contains the Linked List class declaration.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#include "Item.h"

Item::Item() {

}

Item::Item(double cost_, int* cardinality_) {
    this->cost = cost_;
    this->cardinality = cardinality_;
    this->next = NULL;
}

Item::~Item() {
    delete [] this->cardinality;
}
