#include "viterbi.h"
#include <iostream>
#include <math.h>

Viterbi::Viterbi()
{
    //undefined
}

Viterbi::Viterbi(const HMM& hmm)
{
    this->hmm = hmm;
    std::cout << "Hello world from Viterbi!" << std::endl;
}

std::vector<std::string> Viterbi::solve(const std::vector<std::string>& sequence)
{
    std::vector<std::string> state_list = this->hmm.get_state_list();
    std::vector<std::string> obs_list = this->hmm.get_observation_list();
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
        std::string state = state_list[i];
        viterbi[0][i] = log(this->hmm.get_init_prob(state)) + log(this->hmm.get_obs_prob(state, sequence[0]));
        backpointer[0][i] = 0;
    }

    //dynamic programming recursion
    for (int i = 1; i < seq_size; i++)          //for each observation
    {
        for (int j = 0; j < state_size; j++)        //for each current state
        {
            std::string current_state = state_list[j];
            double max_so_far = -9999999;
            int max_so_far_index = 0;
            for (int k = 0; k < state_size; k++)        //go through all previous state
            {
                std::string previous_state = state_list[k];
                double prob = viterbi[i-1][k] + log(this->hmm.get_trans_prob(previous_state, current_state)) + log(this->hmm.get_obs_prob(current_state, sequence[i]));
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
    std::vector<std::string> answer;
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
    std::reverse(answer.begin(), answer.end());

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