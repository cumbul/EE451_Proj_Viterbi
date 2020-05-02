#include <iostream>
#include <mpi.h>
#include "util.h"
#include "viterbi.h"
#include <map>
#include <time.h>

using namespace std;

#define ROOT 0
#define BUFFER_SIZE 20
#define TAG 200

HMM bcast_hmm(HMM& hmm, int my_rank);
vector<string> scatter_seq(HMM& hmm, vector<string>& seq, int my_rank, int num_nodes);
//vector<string> worker_recv_seq();

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

    //distribut data
    hmm = bcast_hmm(hmm, my_rank);

    //***************************debug***************************
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == ROOT)
    {
        cout << "Original sequence:" << endl;
        for (string str : seq)
            cout << str << " ";
        cout << endl;
    }
    //***************************debug***************************
    seq = scatter_seq(hmm, seq, my_rank, num_nodes);

    if (my_rank == ROOT)
    {
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        cout << "Parallel Viterbi takes " << time << " s." << endl;
    }

    //***************************debug***************************
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Node " << my_rank << "sequence:" << endl;
    for (string str : seq)
        cout << str << " ";
    cout << endl;

    cout << "***************************" << endl;
    //***************************debug***************************

    //host_scatter_seq(seq, num_nodes);


    //phase 1
    
    //send last vector to node 1

    //wait until receiving the backtracking result from node 1

    //backtrack


    //workers///

    //once receive necessary data they can start phase 1 right away
    

    //when done phase 1, send the last vector to next node
    
    //loop until converged
        //receive vector from previous node
        //phase 2 : fix vectors
        //if converged, send "converge" signal to root.
        //send last vector to next node

    //if this is the last node, start backtracking, then send the result to previous node.
    //else, receive the backtracking result from next node, backtrack, and send the result to previous node.
    MPI_Finalize();
    
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

//vector<string> worker_recv_seq()
//{

//}
