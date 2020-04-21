#ifndef VITERBI_H
#define VITERBI_H

#include "util.h"

/****************************
Vanilla Viterbi algorithm
****************************/
class Viterbi
{
private:
    HMM hmm;
public:
    Viterbi();
    Viterbi(const HMM& hmm);
    std::vector<std::string> solve(const std::vector<std::string>& sequence);
};

#endif