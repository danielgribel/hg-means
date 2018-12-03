#include "Heap.h"

HeapPdi::HeapPdi() {

}

HeapPdi::~HeapPdi() {
    // ClearHeap();
}

void HeapPdi::PushMax(double cost, int val) {
    heap.push_back( std::pair<double, int>(cost, val) );
    push_heap(heap.begin(), heap.end());
}

int HeapPdi::PopMax() {
    make_heap(heap.begin(), heap.end());
    std::pair<double, int> pair = heap.front();
    pop_heap(heap.begin(), heap.end());
    heap.pop_back();
    return pair.second;
}

std::pair<double, int> HeapPdi::FrontMax() {
    make_heap(heap.begin(), heap.end());
    return heap.front();
}

// void HeapPdi::PushMin(double cost, int val) {
//     heap.push_back( std::pair<double, int>(cost, val) );
//     push_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
// }

// int HeapPdi::PopMin() {
//     make_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
//     std::pair<double, int> pair = heap.front();
//     pop_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
//     heap.pop_back();
//     return pair.second;
// }

// std::pair<double, int> HeapPdi::FrontMin() {
//     make_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
//     return heap.front();
// }

// void HeapPdi::ClearHeap() {
//     std::vector< std::pair<double, int> >().swap(heap);
// }