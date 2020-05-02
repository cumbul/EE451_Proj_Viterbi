#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <vector>
#include <string>

using namespace std;

typedef map<string, map<string, double> > Table;  //a 2d array that uses string as index and it stores double

class HMM
{
private:
    //trans_prob[A][B] = transition probability(double) from previous state A(string) to current state B(string).
    Table trans_prob;

    //obs_prob[A][O] = state observation likelihood(double) of the observation symbol O(string) given the current state A(string).
    Table obs_prob;  

    //init_prob[state] = initial probability of each state
    map<string, double> init_prob;
public:
    HMM();
    HMM(const Table& transition_prob, const Table& observation_prob, const map<string, double>& init_prob);
    double get_trans_prob(const string& previous, const string& next);
    double get_obs_prob(const string& state, const string& observation);
    double get_init_prob(const string& state);
    vector<string> get_state_list();
    vector<string> get_observation_list();
};


namespace Util{
    HMM getToyExample1();
    const vector<string> getToyExample1_Test1();
    const vector<string> getToyExample1_Ans1();
    const vector<string> getToyExample1_Test2();
    const vector<string> getToyExample1_Ans2();

    HMM getToyExample2();
    const vector<string> getToyExample2_Test();
    const vector<string> getToyExample2_Ans();

    HMM getRandomHMM(int num_state, int num_obs);
    vector<string> getRandomSequence(HMM& hmm, int num_steps);
}

#endif
