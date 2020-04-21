#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <vector>

typedef std::map<std::string, std::map<std::string, double> > Table;  //a 2d array that uses string as index and it stores double

class HMM
{
private:
    //trans_prob[A][B] = transition probability(double) from previous state A(string) to current state B(string).
    Table trans_prob;

    //obs_prob[A][O] = state observation likelihood(double) of the observation symbol O(string) given the current state A(string).
    Table obs_prob;  

    //init_prob[state] = initial probability of each state
    std::map<std::string, double> init_prob;
public:
    HMM();
    HMM(const Table& transition_prob, const Table& observation_prob, const std::map<std::string, double>& init_prob);
    double get_trans_prob(const std::string& previous, const std::string& next);
    double get_obs_prob(const std::string& state, const std::string& observation);
    double get_init_prob(const std::string& state);
    std::vector<std::string> get_state_list();
    std::vector<std::string> get_observation_list();
};


namespace Util{
    HMM getToyExample1();
    const std::vector<std::string> getToyExample1_Test1();
    const std::vector<std::string> getToyExample1_Ans1();
    const std::vector<std::string> getToyExample1_Test2();
    const std::vector<std::string> getToyExample1_Ans2();
}

#endif