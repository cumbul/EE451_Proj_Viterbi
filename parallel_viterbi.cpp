#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cstring>
#include <omp.h>
#include <algorithm>
#include "parallel_viterbi.h"

/****************************
Parallel Viterbi Algorithm using Linear Tropical Dynamic Programming
****************************/
LTDPViterbi::LTDPViterbi(const HMM& hmm, int num_processor) : ParallelViterbi(hmm, num_processor)
{
    srand(time(NULL));
    cout << "LTDP Viterbi: Using " << num_processor << " processors." << endl;
}

vector<string> LTDPViterbi::solve(const vector<string>& sequence)
{
    vector<string> state_list = hmm.get_state_list();
    vector<string> obs_list = hmm.get_observation_list();
    int seq_size = sequence.size();
    int state_size = state_list.size();
    bool conv[num_processor];

    //initialization
    double** viterbi = new double*[seq_size];
    int** pred = new int*[seq_size];
    for (int i = 0; i < seq_size; i++)
    {
        viterbi[i] = new double[state_size];
        pred[i] = new int[state_size];
    }

    for (int i = 0; i < seq_size; i++)
    	for (int j = 0; j < state_size; j++)
	        pred[i][j] = 0;

    for (int i = 0; i < num_processor; i++) conv[i] = false;

    //parallel forward phase
    #pragma omp parallel num_threads(num_processor) shared(state_list, obs_list, seq_size, state_size, viterbi, pred)
    {
        int tid = omp_get_thread_num()+1;
        double previous_stage[state_size], hold[state_size];
        int left_p = (seq_size / num_processor) * (tid - 1);
        int right_p = min((seq_size / num_processor) * tid, seq_size-1);
	
        //Only the first processor gets the true initial probabilities. Other processors will get a random vector.
        if (tid == 1)
        {
            for (int i = 0; i < state_size; i++)
            {
                string state = state_list[i];
                viterbi[0][i] = previous_stage[i] = log(hmm.get_init_prob(state)) + log(hmm.get_obs_prob(state, sequence[0]));
                pred[0][i] = 0;
            }
        }
        else
        {
            for (int i = 0; i < state_size; i++)
                previous_stage[i] = (double)(rand() % 100) / 100.0;
        }
        
        for (int i = left_p + 1; i <= right_p; i++)
        {
            for (int j = 0; j < state_size; j++)        //for each current state
            {
                string current_state = state_list[j];
                double max_so_far = -INFINITY;
                int max_so_far_index = 0;
                for (int k = 0; k < state_size; k++)        //go through all previous state
                {
                    string previous_state = state_list[k];
                    double prob = previous_stage[k] + log(hmm.get_trans_prob(previous_state, current_state)) + log(hmm.get_obs_prob(current_state, sequence[i]));
                    if (prob > max_so_far)
                    {
                        max_so_far = prob;
                        max_so_far_index = k;
                    }
                }
                viterbi[i][j] = hold[j] = max_so_far;
                pred[i][j] = max_so_far_index;
            }
            memcpy(previous_stage, hold, sizeof(double) * state_size);
        }
    }

    if (num_processor > 1)
    {
        // till  convergence (fix up loop)
        bool converged = false;
        do
        {
            #pragma omp parallel num_threads(num_processor-1) shared(state_list, obs_list, seq_size, state_size, viterbi, pred, conv)
            {
                int tid = omp_get_thread_num() + 2;
                conv[tid-2] = false;
                int left_p = (seq_size / num_processor) * (tid - 1);
                int right_p = std::min((seq_size / num_processor) * tid, seq_size-1);

                //obtain final solution from previous processor
                double s[state_size], hold[state_size];
                memcpy(s, viterbi[left_p], sizeof(double) * state_size);

                for (int i = left_p + 1; i <= right_p; i++)
                {
                    for (int j = 0; j < state_size; j++)        //for each current state
                    {
                        string current_state = state_list[j];
                        double max_so_far = -INFINITY;
                        int max_so_far_index = 0;
                        for (int k = 0; k < state_size; k++)        //go through all previous state
                        {
                            string previous_state = state_list[k];
                            double prob = s[k] + log(hmm.get_trans_prob(previous_state, current_state)) + log(hmm.get_obs_prob(current_state, sequence[i]));
                            if (prob > max_so_far)
                            {
                                max_so_far = prob;
                                max_so_far_index = k;
                            }
                        }
			            hold[j] = max_so_far;
                        pred[i][j] = max_so_far_index;
                    }

		            memcpy(s, hold, sizeof(double) * state_size);

                    if (_isParallel(s, viterbi[i], state_size))
                    {
                        conv[tid-2] = true;
                        break;
                    }
                }
            }
            converged = conv[0];
            for (int i = 0; i < num_processor-1; i++)
                converged = converged && conv[i];
        } while (!converged);
    }

    //backtracking
    vector<string> answer(seq_size);
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
        answer[i-1] = state_list[bestpathpointer];
        bestpathpointer = pred[i-1][bestpathpointer];
    }

    //clean up
    for (int i = 0; i < seq_size; i++)
    {
        delete [] viterbi[i];
        delete [] pred[i];
    }
    delete [] viterbi;
    delete [] pred;

    return answer;
}

bool LTDPViterbi::_isParallel(double* A, double* B, int size)
{
    double diff = A[0] - B[0];
    double delta;
    for (int i = 1; i < size; i++)
    {
        delta = A[i] - B[i];
        if (fabs(diff - delta) > EPSILON)
            return false;
    }
    return true;
}

/****************************
Parallel Viterbi Algorithm using LTDP with "Processor Oblivious"
****************************/
vector<string> LTDPViterbi_Oblivious::solve(const vector<string>& sequence)
{
    //if num_processor = 1, then it is same as LTDPViterbi
    if (num_processor == 1) return LTDPViterbi::solve(sequence);

    vector<string> state_list = hmm.get_state_list();
    vector<string> obs_list = hmm.get_observation_list();
    int seq_size = sequence.size();
    int state_size = state_list.size();
    int t = seq_size;
    int c = pow(2, max(int(log2(t)) - 4, 2));
    int size = c;
    int q = t / size;

    //initialization
    double** viterbi = new double*[seq_size];
    int** pred = new int*[seq_size];
    for (int i = 0; i < seq_size; i++)
    {
        viterbi[i] = new double[state_size];
        pred[i] = new int[state_size];
    }

    for (int i = 0; i < seq_size; i++)
    	for (int j = 0; j < state_size; j++)
	        pred[i][j] = 0;

    //forward phase
    #pragma omp parallel num_threads(q) shared(state_list, obs_list, seq_size, state_size, viterbi, pred)
    {
        int tid = omp_get_thread_num();
        double previous_stage[state_size], hold[state_size];
        int left_p  = tid * size;
        int right_p = min(left_p + size, seq_size - 1);

        //Only the first processor gets the true initial probabilities. Other processors will get a random vector. 
        if (tid == 0)
        {
            for (int i = 0; i < state_size; i++)
            {
                string state = state_list[i];
                viterbi[0][i] = previous_stage[i] = log(hmm.get_init_prob(state)) + log(hmm.get_obs_prob(state, sequence[0]));
                pred[0][i] = 0;
            }
        }
        else
        {
            for (int i = 0; i < state_size; i++)
                previous_stage[i] = (double)(rand() % 100) / 100.0;
        }

        for (int i = left_p + 1; i <= right_p; i++)
        {
            for (int j = 0; j < state_size; j++)        //for each current state
            {
                string current_state = state_list[j];
                double max_so_far = -INFINITY;
                int max_so_far_index = 0;
                for (int k = 0; k < state_size; k++)        //go through all previous state
                {
                    string previous_state = state_list[k];
                    double prob = previous_stage[k] + log(hmm.get_trans_prob(previous_state, current_state)) + log(hmm.get_obs_prob(current_state, sequence[i]));
                    if (prob > max_so_far)
                    {
                        max_so_far = prob;
                        max_so_far_index = k;
                    }
                }
                viterbi[i][j] = hold[j] = max_so_far;
                pred[i][j] = max_so_far_index;
            }
            memcpy(previous_stage, hold, sizeof(double) * state_size);            
        }
    }
    //fixup phase
    bool converged = false;
    for (int bound = int(log2(c)); bound < int(log2(t)) - 1 && !converged; bound++)
    {
        size = pow(2, bound);
        q = t / (2 * size);

        bool* conv = new bool[q];
        for (int i = 0; i < q; i++) conv[i] = false;

        #pragma omp parallel num_threads(q) shared(state_list, obs_list, seq_size, state_size, viterbi, pred, conv)
        {
            int tid = omp_get_thread_num();
            int left_p = (((2 * tid) + 1) * size);
            int right_p = min(left_p, seq_size - 1);
            double u[state_size], hold[state_size];
            memcpy(u, viterbi[left_p], sizeof(double) * state_size);

            for (int i = left_p + 1; i <= right_p; i++)
            {
                for (int j = 0; j < state_size; j++)        //for each current state
                {
                    string current_state = state_list[j];
                    double max_so_far = -INFINITY;
                    int max_so_far_index = 0;
                    for (int k = 0; k < state_size; k++)        //go through all previous state
                    {
                        string previous_state = state_list[k];
                        double prob = u[k] + log(hmm.get_trans_prob(previous_state, current_state)) + log(hmm.get_obs_prob(current_state, sequence[i]));
                        if (prob > max_so_far)
                        {
                            max_so_far = prob;
                            max_so_far_index = k;
                        }
                    }
			        hold[j] = max_so_far;
                    pred[i][j] = max_so_far_index;
                }
                memcpy(u, hold, sizeof(double) * state_size);
                if (_isParallel(u, viterbi[i], state_size))
                {
                    conv[tid] = true;
                    break;
                }
		memcpy(viterbi[i], hold, sizeof(double) * state_size);    
            }
        }
        converged = conv[0];
        for (int i = 0; i < q; i++) converged = conv[i] && converged;
        delete [] conv;
    }

    //backtracking
    vector<string> answer(seq_size);
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
        answer[i-1] = state_list[bestpathpointer];
        bestpathpointer = pred[i-1][bestpathpointer];
    }

    //clean up
    for (int i = 0; i < seq_size; i++)
    {
        delete [] viterbi[i];
        delete [] pred[i];
    }
    delete [] viterbi;
    delete [] pred;

    return answer;
}
