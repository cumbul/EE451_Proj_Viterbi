#include <iostream>
#include <mpi.h>
#include "util.h"
#include "viterbi.h"
#include <map>
#include <time.h>

using namespace std;

#define ROOT 0
#define BUFFER_SIZE 30
#define TAG 200

void host_distribute_hmm(HMM& hmm, int num_nodes);
//void host_scatter_seq(vector<string>& seq, int num_nodes);
HMM worker_recv_hmm();
//vector<string> worker_recv_seq();

int main(int argc, char *argv[])
{
    //change the parameter here
    int num_state = 4, num_obs = 4, seq_length = 4;

    MPI_Init(&argc, &argv);
    int num_nodes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == ROOT)
    {
        HMM hmm = Util::getRandomHMM(num_state, num_obs);
        vector<string> seq = Util::getRandomSequence(hmm, seq_length);
        struct timespec start, stop;
        double time;

        //***************************debug***************************
        cout << "Debug Host:" << endl;
        vector<string> state_list = hmm.get_state_list();
        cout << "state trans" << endl;
        for (int i = 0; i < state_list.size(); i++)
        {
            for (int j = 0; j < state_list.size(); j++)
            {
                cout << state_list[i] << " -> " << state_list[j] << " : " << hmm.get_trans_prob(state_list[i], state_list[j]) << endl;
            }
        }
        cout << endl;
        cout << "obs emission" << endl;
        vector<string> obs_list = hmm.get_observation_list();
        for (int i = 0; i < state_list.size(); i++)
        {
            for (int j = 0; j < obs_list.size(); j++)
            {
                cout << state_list[i] << " -> " << obs_list[j] << " : " << hmm.get_obs_prob(state_list[i], obs_list[j]) << endl;
            }
        }
        cout << "***************************" << endl;
        //***************************debug***************************


        //Vanilla Viterbi is only computed in node 0
        cout << "Vanilla Viterbi: " << endl;
        Viterbi viterbi(hmm);
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
        vector<string> vanilla_ans = viterbi.solve(seq);
        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        cout << "Vanilla Viterbi takes " << time << "s." << endl;  

        cout << "Parallel Viterbi: " << endl;
        if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}

        //distribut data
        host_distribute_hmm(hmm, num_nodes);
        //host_scatter_seq(seq, num_nodes);


        //phase 1
        
        //send last vector to node 1

        //wait until receiving the backtracking result from node 1

        //backtrack

        if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
        time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
        cout << "Parallel Viterbi takes " << time << " s." << endl;
    }
    else
    {
        //receive data
        cout << "hello world from worker " << my_rank << "!" << endl; 
        HMM hmm = worker_recv_hmm();
        //vector<string> seq = worker_recv_seq();

        //***************************debug***************************
        cout << "Debug Worker " << my_rank << ":" << endl;
        vector<string> state_list = hmm.get_state_list();
        cout << "Worker " << my_rank << " state trans:" << endl;
        for (int i = 0; i < state_list.size(); i++)
        {
            for (int j = 0; j < state_list.size(); j++)
            {
                cout << state_list[i] << " -> " << state_list[j] << " : " << hmm.get_trans_prob(state_list[i], state_list[j]) << endl;
            }
        }
        cout << endl;
        cout << "Worker " << my_rank << " obs emission:" << endl;
        vector<string> obs_list = hmm.get_observation_list();
        for (int i = 0; i < state_list.size(); i++)
        {
            for (int j = 0; j < obs_list.size(); j++)
            {
                cout << state_list[i] << " -> " << obs_list[j] << " : " << hmm.get_obs_prob(state_list[i], obs_list[j]) << endl;
            }
        }
        cout << "***************************" << endl;
        //***************************debug***************************

        //once receive necessary data they can start phase 1 right away
        

        //when done phase 1, send the last vector to next node
        
        //loop until converged
            //receive vector from previous node
            //phase 2 : fix vectors
            //if converged, send "converge" signal to root.
            //send last vector to next node

        //if this is the last node, start backtracking, then send the result to previous node.
        //else, receive the backtracking result from next node, backtrack, and send the result to previous node.
    }
    MPI_Finalize();
    
    return 0;
}

void host_distribute_hmm(HMM& hmm, int num_nodes)
{
    vector<string> state_list = hmm.get_state_list();
    //send number of states
    int num_states = state_list.size();
    MPI_Bcast(&num_states, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	cout << "Host: sent num_states!" << endl;
    //send state_list
    for (int i = 0; i < num_states; i++)
    {
        MPI_Bcast(&state_list[i][0], state_list[i].size(), MPI_CHAR, ROOT, MPI_COMM_WORLD);
    }

    //send state transition probabilities
    for (int i = 0; i < num_states; i++)
    {
        string previous = state_list[i];
        for (int j = 0; j < num_states; j++)
        {
            string next = state_list[j];
            double prob = hmm.get_trans_prob(previous, next);
            MPI_Bcast(&prob, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        }
    }

    //send number of observations
    vector<string> obs_list = hmm.get_observation_list();
    int num_obs = obs_list.size();
    MPI_Bcast(&num_obs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    //send observation list
    for (int i = 0; i < num_obs; i++)
    {
        MPI_Bcast(&obs_list[i][0], obs_list[i].size(), MPI_CHAR, ROOT, MPI_COMM_WORLD);
    }

    //send observation probabilities
    for (int i = 0; i < num_states; i++)
    {
        string state = state_list[i];
        for (int j = 0; j < num_obs; j++)
        {
            string obs = obs_list[j];
            double prob = hmm.get_obs_prob(state, obs);
            MPI_Bcast(&prob, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        }
    }
}

//void host_scatter_seq(vector<string>& seq, int num_nodes)
//{
//}

HMM worker_recv_hmm()
{
	cout << "waiting for num_states..." << endl;
    Table trans_prob, obs_prob;
    //receive number of states
    int num_states;
    MPI_Recv(&num_states, 1, MPI_INT, ROOT, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    vector<string> state_list(num_states);
	cout << "received num_state!" << endl;
    //receive state_list
    char buffer[BUFFER_SIZE];
    for (int i = 0; i < num_states; i++)
    {
        MPI_Bcast(&buffer, BUFFER_SIZE, MPI_CHAR, ROOT, MPI_COMM_WORLD);
        state_list[i] = string(buffer);
    }
	cout << "done state_list!" << endl;
    //receive state transition probabilities
    //send state transition probabilities
    for (int i = 0; i < num_states; i++)
    {
        string previous = state_list[i];
        for (int j = 0; j < num_states; j++)
        {
            string next = state_list[j];
            double prob;
            MPI_Bcast(&prob, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            trans_prob[previous][next] = prob;
        }
    }
	cout << "done trans_prob!" << endl;
    //receive number of observation
    int num_obs;
    MPI_Bcast(&num_obs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    vector<string> obs_list(num_obs);

    //receive observation list
    for (int i = 0; i < num_obs; i++)
    {
        MPI_Bcast(&buffer, BUFFER_SIZE, MPI_CHAR, ROOT, MPI_COMM_WORLD);
        obs_list[i] = string(buffer);
    }
	cout << "done obs_list!" << endl;
    //receive observation probabilities
    for (int i = 0; i < num_states; i++)
    {
        string state = state_list[i];
        for (int j = 0; j < num_obs; j++)
        {
            string obs = obs_list[j];
            double prob;
            MPI_Bcast(&prob, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            obs_prob[state][obs] = prob;
        }
    }
	cout << "done obs probs!" << endl;
    return HMM(trans_prob, obs_prob, map<string, double>());
}

//vector<string> worker_recv_seq()
//{

//}
