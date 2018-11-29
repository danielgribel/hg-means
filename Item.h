#ifndef Item_h
#define Item_h

#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

class Item {

    public:

        double cost;
        
        int* cardinality;
        
        Item * next;

        Item();
        
        Item(double cost_, int* cardinality_);
        
        ~Item();
};

#endif