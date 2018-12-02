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
    
        void PushMax(double cost, int val);
    
        int PopMax();
    
        void PushMin(double cost, int val);
    
        int PopMin();

        std::pair<double, int> FrontMax();
	
    	std::pair<double, int> FrontMin();

        std::vector< std::pair<double, int> > GetHeap() { return heap; };
    
        void ClearHeap();
};

#endif