#include <iostream>

#include "viterbi.h"
#include "parallel_viterbi.h"
#include "util.h"
#include <map>
#include <time.h>

using namespace std;
//typedef map<string, map<string, float> > Table; 

int main()
{
    //change the parameter here
    int num_state = 40, num_obs = 80, seq_length = 1024, num_cores = 2;

    HMM hmm = Util::getRandomHMM(num_state, num_obs);
    vector<string> seq = Util::getRandomSequence(hmm, seq_length);
    struct timespec start, stop;
    double time;

    cout << "Vanilla Viterbi: " << endl;

    Viterbi viterbi(hmm);
    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    vector<string> vanilla_ans = viterbi.solve(seq);
    if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
    time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
    cout << "Vanilla Viterbi takes " << time << "s." << endl;    
    
    cout << "Parallel Viterbi: " << endl;
    LTDPViterbi_Oblivious ltdp_viterbi(hmm, num_cores);
    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) { perror( "clock gettime" );}
    vector<string> ltdp_ans = ltdp_viterbi.solve(seq);
    if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) { perror( "clock gettime" );}
    time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;
    cout << "Parallel Viterbi takes " << time << " s." << endl;

    cout << "Checking correctness...";
    bool correct = true;
    for (int i = 0; i < seq_length; i++)
    {
	if (vanilla_ans[i].compare(ltdp_ans[i]) != 0)
	{
	    correct = false;
	    break;
	}
    }
    if (correct)
	cout << "Correct!";
    else
	cout << "Incorrect!";
    cout << endl;
}
