#include "PbData.h"

PbData::PbData() {

}

PbData::PbData(double* data, int n, int d, int m) {
    this->data = data;
    this->n = n;
    this->d = d;
    this->m = m;
}

PbData::~PbData() {

}