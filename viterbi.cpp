#include "viterbi.h"
#include <iostream>
#include <math.h>
#include <algorithm>

Viterbi::Viterbi()
{
    //undefined
}

Viterbi::Viterbi(const HMM& hmm)
{
    this->hmm = hmm;
}

vector<string> Viterbi::solve(const vector<string>& sequence)
{
    vector<string> state_list = hmm.get_state_list();
    vector<string> obs_list = hmm.get_observation_list();
    int seq_size = sequence.size();
    int state_size = state_list.size();

    //initialization
    double** viterbi = new double*[seq_size];
    int** backpointer = new int*[seq_size];
    for (int i = 0; i < seq_size; i++)
    {
        viterbi[i] = new double[state_size];
        backpointer[i] = new int[state_size];
    }

    for (int i = 0; i < state_size; i++)
    {
        string state = state_list[i];
        viterbi[0][i] = log(hmm.get_init_prob(state)) + log(hmm.get_obs_prob(state, sequence[0]));
        backpointer[0][i] = 0;
    }

    //dynamic programming recursion
    for (int i = 1; i < seq_size; i++)          //for each observation
    {
        for (int j = 0; j < state_size; j++)        //for each current state
        {
            string current_state = state_list[j];
            double max_so_far = -INFINITY;
            int max_so_far_index = 0;
            for (int k = 0; k < state_size; k++)        //go through all previous state
            {
                string previous_state = state_list[k];
                double prob = viterbi[i-1][k] + log(hmm.get_trans_prob(previous_state, current_state)) + log(hmm.get_obs_prob(current_state, sequence[i]));
                if (prob > max_so_far)
                {
                    max_so_far = prob;
                    max_so_far_index = k;
                }
            }
            viterbi[i][j] = max_so_far;
            backpointer[i][j] = max_so_far_index;
        }
    }

    //backtracking
    vector<string> answer;
    double bestprob = viterbi[seq_size-1][0];
    int bestpathpointer = 0;
    for (int i = 0; i < state_size; i++)
    {
        if (viterbi[seq_size-1][i] > bestprob)
        {
            bestprob = viterbi[seq_size-1][i];
            bestpathpointer = i;
        }
    }
    
    for (int i = seq_size; i > 0; i--)
    {
        answer.push_back(state_list[bestpathpointer]);
        bestpathpointer = backpointer[i-1][bestpathpointer];
    }
    reverse(answer.begin(), answer.end());

    //clean up
    for (int i = 0; i < seq_size; i++)
    {
        delete [] viterbi[i];
        delete [] backpointer[i];
    }
    delete [] viterbi;
    delete [] backpointer;

    return answer;
}