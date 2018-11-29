#ifndef Heap_Pdi_H
#define Heap_Pdi_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <functional>

class HeapPdi {

    private:
    
        std::vector< std::pair<double, int> > heap;
    
    public:
    
        HeapPdi();
    
        ~HeapPdi();
    
        void push_max(double cost, int val);
    
        int pop_max();
    
        void push_min(double cost, int val);
    
        int pop_min();
    
        std::vector< std::pair<double, int> > getHeap() { return heap; };
    
        void setHeap(std::vector< std::pair<double, int> > aHeap);
    
        void clearHeap();
    
        std::pair<double, int> front_max();
	
    	std::pair<double, int> front_min();
};

#endif