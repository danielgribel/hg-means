#include "PbRun.h"

PbRun::PbRun() {
    this->solution = NULL;
    this->time = 0.0;
    this->lastImprovement = 0;
    this->nbIter = 0.0;
    this->avgAlpha = 0.0;
    this->isValid = false;
}

PbRun::PbRun(Solution* solution, double time, int lastImprovement, double nbIter, double avgAlpha, bool isValid) {
    this->solution = solution;
    this->time = time;
    this->lastImprovement = lastImprovement;
    this->nbIter = nbIter;
    this->avgAlpha = avgAlpha;
    this->isValid = isValid;
}

PbRun::~PbRun() {

}