#include <iostream>
#include "parallel_viterbi.h"

/****************************
Parallel Viterbi Algorithm using Linear Tropical Dynamic Programming
****************************/
LTDPViterbi::LTDPViterbi(const HMM& hmm, int num_processor) : ParallelViterbi(hmm, num_processor)
{
    cout << "Hello world from LTDP Viterbi!" << endl;
    cout << "Using " << num_processor << " processors." << endl;
}

vector<string> LTDPViterbi::solve(const vector<string>& sequence)
{
    //vector<string> state_list = hmm.get_state_list();
    //vector<string> obs_list = hmm.get_observation_list();
    //int seq_size = sequence.size();
    //int state_size = state_list.size();
    vector<string> hi;
    hi.push_back("not implemented yet");
    return hi;
}