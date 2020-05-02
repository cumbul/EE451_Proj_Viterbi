#include <iostream>
#include <map>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "util.h"
#include "viterbi.h"

using namespace std;

#define ROOT 0
#define BUFFER_SIZE 20
#define TAG 200

HMM bcast_hmm(HMM& hmm, int my_rank);
vector<string> scatter_seq(HMM& hmm, vector<string>& seq, int my_rank, int num_nodes);

int main(int argc, char *argv[])
{
    //change the parameter here
    int num_state = 4, num_obs = 4, seq_length = 16;

    MPI_Init(&argc, &argv);
    int num_nodes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (seq_length % num_nodes != 0)
    {
        cout << "Error: The seq_length must be divisible to number of nodes!" << endl;
        MPI_Finalize();
        return 0;
    }

    HMM hmm;
    vector<string> seq, vanilla_ans;
    struct timespec start, stop;
    double time;

    //******************************vanilla******************************
    if (my_rank == ROOT)
    {
        hmm = Util::getRandomHMM(num_state, num_obs);
        seq = Util::getRandomSequence(hmm, seq_length);

        //Vanilla Viterbi is only computed in node 0
        cout << "Vanilla Viterbi: " << endl;
        Viterbi viterbi(hmm);
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
        vanilla_ans = viterbi.solve(seq);
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        cout << "Vanilla Viterbi takes " << time << "s." << endl;  

        cout << "Parallel Viterbi: " << endl;
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    }
    //******************************vanilla******************************

    //******************************distribute data******************************
    hmm = bcast_hmm(hmm, my_rank);
    seq = scatter_seq(hmm, seq, my_rank, num_nodes);
    //******************************distribute data******************************

    //******************************forward phase******************************
    vector<string> state_list = hmm.get_state_list();
    vector<string> obs_list = hmm.get_observation_list();
    int seq_size = seq.size();
    int state_size = state_list.size();
    double** viterbi = new double*[seq_size];
    int** pred = new int*[seq_size];
    double previous_stage[state_size], hold[state_size];

    //initialization
    for (int i = 0; i < seq_size; i++)
    {
        viterbi[i] = new double[state_size];
        pred[i] = new int[state_size];
    }

    for (int i = 0; i < seq_size; i++)
    	for (int j = 0; j < state_size; j++)
	        pred[i][j] = 0;

    //Only the first node gets the true initial probabilities. Other node will get a random vector.
    if (my_rank == ROOT)
    {
        for (int i = 0; i < state_size; i++)
        {
            string state = state_list[i];
            viterbi[0][i] = previous_stage[i] = log(hmm.get_init_prob(state)) + log(hmm.get_obs_prob(state, seq[0]));
            pred[0][i] = 0;
        }
    }
    else
    {
        for (int i = 0; i < state_size; i++)
            previous_stage[i] = (double)(rand() % 100) / 100.0;
    }

    for (int i = 1; i < seq.size(); i++)
    {
        for (int j = 0; j < state_size; j++)
        {
            string current_state = state_list[j];
            double max_so_far = -INFINITY;
            int max_so_far_index = 0; 
            for (int k = 0; k < state_size; k++)
            {
                string previous_state = state_list[k];
                double prob = previous_stage[k] + log(hmm.get_trans_prob(previous_state, current_state)) + log(hmm.get_obs_prob(current_state, seq[i]));
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

    //******************************forward phase******************************

    //******************************fixing phase******************************
    //******************************fixing phase******************************

    //******************************backtrack phase******************************
    //******************************backtrack phase******************************

    //runtime output
    if (my_rank == ROOT)
    {
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        cout << "Parallel Viterbi takes " << time << " s." << endl;
    }

    MPI_Finalize();

    //clean up
    for (int i = 0; i < seq_size; i++)
    {
        delete [] viterbi[i];
        delete [] pred[i];
    }
    delete [] viterbi;
    delete [] pred;
    
    return 0;
}

HMM bcast_hmm(HMM& hmm, int my_rank)
{
    Table trans_prob, obs_prob;
    vector<string> state_list, obs_list;
    int num_states, num_obs;
    char buffer[BUFFER_SIZE];

    if (my_rank == ROOT)
    {
        state_list = hmm.get_state_list();
        obs_list = hmm.get_observation_list();
        num_states = state_list.size();
        num_obs = obs_list.size();
    }

    //send number of states
    MPI_Bcast(&num_states, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    //send state_list
    for (int i = 0; i < num_states; i++)
    {
        if (my_rank == ROOT) strcpy(buffer, &state_list[i][0]);
        MPI_Bcast(buffer, BUFFER_SIZE, MPI_CHAR, ROOT, MPI_COMM_WORLD);
        if (my_rank != ROOT) state_list.push_back(string(buffer));
    }
    //send state transition probabilities
    for (int i = 0; i < num_states; i++)
    {
        string previous = state_list[i];
        for (int j = 0; j < num_states; j++)
        {
            string next = state_list[j];
            double prob = (my_rank == ROOT)? hmm.get_trans_prob(previous, next): 0;
            MPI_Bcast(&prob, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            trans_prob[previous][next] = prob;
        }
    }

    //send number of observations
    MPI_Bcast(&num_obs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    //send observation list
    for (int i = 0; i < num_obs; i++)
    {
        if (my_rank == ROOT) strcpy(buffer, &obs_list[i][0]);
        MPI_Bcast(buffer, BUFFER_SIZE, MPI_CHAR, ROOT, MPI_COMM_WORLD);
        if (my_rank != ROOT) obs_list.push_back(string(buffer));
    }
    //send observation probabilities
    for (int i = 0; i < num_states; i++)
    {
        string state = state_list[i];
        for (int j = 0; j < num_obs; j++)
        {
            string obs = obs_list[j];
            double prob = (my_rank == ROOT)? hmm.get_obs_prob(state, obs): 0;
            MPI_Bcast(&prob, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            obs_prob[state][obs] = prob;
        }
    }    
        
    if (my_rank == 0)
        return hmm;
    else
        return HMM(trans_prob, obs_prob, map<string, double>());
}

vector<string> scatter_seq(HMM& hmm, vector<string>& seq, int my_rank, int num_nodes)
{
    int seq_length = seq.size();
    MPI_Bcast(&seq_length, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    map<string, int> obs_to_int;
    const vector<string>& obs_list = hmm.get_observation_list();
    int element_per_proc = seq_length / num_nodes;

    int* obs_sent = new int[seq.size()];
    int* obs_recv = new int[element_per_proc];
    if (my_rank == ROOT)
    {
        for (int i = 0; i < obs_list.size(); i++)
            obs_to_int[obs_list[i]] = i;
        for (int i = 0; i < seq.size(); i++) 
            obs_sent[i] = obs_to_int[seq[i]];
    }
    MPI_Scatter(obs_sent, element_per_proc, MPI_INT, obs_recv, element_per_proc, MPI_INT, ROOT, MPI_COMM_WORLD);
    vector<string> result(element_per_proc);
    for (int i = 0; i < element_per_proc; i++)
        result[i] = obs_list[obs_recv[i]];
    delete [] obs_sent;
    delete [] obs_recv;
    return result;
}
