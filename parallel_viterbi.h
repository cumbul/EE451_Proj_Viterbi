#ifndef PARALLEL_VITERBI_H
#define PARALLEL_VITERBI_H

#include "viterbi.h"
#include "util.h"

/****************************
Abstract Parallel Viterbi algorithm
****************************/
class ParallelViterbi : public Viterbi
{
protected:
    int num_processor;
public:
    ParallelViterbi() : num_processor(1) {};
    ParallelViterbi(const HMM& hmm) : Viterbi(hmm), num_processor(1) {};
    ParallelViterbi(const HMM& hmm, int num_processor) : Viterbi(hmm), num_processor(num_processor) {};
};

/****************************
Parallel Viterbi Algorithm with Linear Tropical Dynamic Programming
****************************/
class LTDPViterbi : public ParallelViterbi
{
protected:
    const double EPSILON = 0.0001;     //error allowed during calculating parallel vectors
    bool _isParallel(double* A, double* B, int size);
public:
    LTDPViterbi() : ParallelViterbi() {};
    LTDPViterbi(const HMM& hmm) : ParallelViterbi(hmm) {};
    LTDPViterbi(const HMM& hmm, int num_processor);
    virtual vector<string> solve(const vector<string>& sequence) override;
};

/****************************
Parallel Viterbi Algorithm using LTDP with "Processor Oblivious"
****************************/
class LTDPViterbi_Oblivious : public LTDPViterbi
{
public:
    LTDPViterbi_Oblivious() : LTDPViterbi() {};
    LTDPViterbi_Oblivious(const HMM& hmm) : LTDPViterbi(hmm) {};
    LTDPViterbi_Oblivious(const HMM& hmm, int num_processor) : LTDPViterbi(hmm, num_processor) {};
    virtual vector<string> solve(const vector<string>& sequence) override;
};
#endif
