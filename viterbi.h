#ifndef VITERBI_H
#define VITERBI_H

#include "util.h"

/****************************
Vanilla Viterbi algorithm
****************************/
class Viterbi
{
protected:
    HMM hmm;
public:
    Viterbi();
    Viterbi(const HMM& hmm);
    virtual vector<string> solve(const vector<string>& sequence);
};

#endif