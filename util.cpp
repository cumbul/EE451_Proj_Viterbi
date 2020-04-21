#include "util.h"
#include <iostream>
#include <string>
/****************************
Hidden Markov Model 
****************************/
HMM::HMM()
{
    //undefined
}

HMM::HMM(const Table& transition_prob, const Table& observation_prob, const std::map<std::string, double>& init_prob)
{
    this->trans_prob = transition_prob;
    this->obs_prob = observation_prob;
    this->init_prob = init_prob;
}

/*
Get transition probability from previous state to the next state.
*/
double HMM::get_trans_prob(const std::string& previous, const std::string& next)
{
    //If the input not found in table, show warning. To be deleted after debugging.
    if (this->trans_prob.find(previous) == this->trans_prob.end())
        std::cout << "Warning: previous argument '" << previous << "' in HMM::get_trans_prob() does not exist. Return 0 by default." << std::endl; 
    else if (this->trans_prob[previous].find(next) == this->trans_prob[previous].end())
        std::cout << "Warning: next argument '" << next << "' in HMM::get_trans_prob() does not exist. Return 0 by default." << std::endl; 

    return this->trans_prob[previous][next];
}

/*
Get observation likelihood of the observation given the current state.
*/
double HMM::get_obs_prob(const std::string& state, const std::string& observation)
{
    //If the input not found in table, show warning. To be deleted after debugging.
    if (this->obs_prob.find(state) == this->obs_prob.end())
        std::cout << "Warning: state argument '" << state << "' in HMM::get_obs_prob() does not exist. Return 0 by default." << std::endl; 
    else if (this->obs_prob[state].find(observation) == this->obs_prob[state].end())
        std::cout << "Warning: observation argument '" << observation << "' in HMM::get_obs_prob() does not exist. Return 0 by default." << std::endl; 

    return this->obs_prob[state][observation];
}

double HMM::get_init_prob(const std::string& state)
{
    //If the input not found in table, show warning. To be deleted after debugging.
    if (this->init_prob.find(state) == this->init_prob.end())
        std::cout << "Warning: state argument '" << state << "' in HMM::get_init_prob() does not exist. Return 0 by default." << std::endl;
    return this->init_prob[state];
}

/*
Get a list of all possible states in HMM.
*/
std::vector<std::string> HMM::get_state_list()
{
    std::vector<std::string> result;
    for(Table::iterator it = this->trans_prob.begin(); it != this->trans_prob.end(); it++)
    {
        result.push_back(it->first);
    }
    return result;
}

/*
Get a list of all possible observation in HMM.
*/
std::vector<std::string> HMM::get_observation_list()
{
    std::vector<std::string> result;
    std::map<std::string, double> first_table = this->obs_prob.begin()->second;
    for(std::map<std::string, double>::iterator it = first_table.begin(); it != first_table.end(); it++)
    {
        result.push_back(it->first);
    }
    return result;
}

/****************************
Toy HMM model for testing
****************************/
/* return a toy HMM model as shown in www.cis.upenn.edu/~cis262/notes/Example-Viterbi-DNA.pdf */
HMM Util::getToyExample1()
{
    Table test_trans_prob{
		{"H",
            {
                {"H", 0.5},
                {"L", 0.5}
            }
        },
		{"L",
            {
                {"H", 0.4},
                {"L", 0.6}
            }
        }
    };   
    Table test_obs_prob{
        {"H",
            {
                {"A", 0.2},
                {"C", 0.3},
                {"G", 0.3},
                {"T", 0.2}
            }
        },
        {"L",
            {
                {"A", 0.3},
                {"C", 0.2},
                {"G", 0.2},
                {"T", 0.3}
            }
        }
    };
    std::map<std::string, double> init_prob{
        {"H", 0.5},
        {"L", 0.5}
    };
    return HMM(test_trans_prob, test_obs_prob, init_prob);
}

const std::vector<std::string> Util::getToyExample1_Test1()
{
    std::string test = "GGCACTGAA";
    std::vector<std::string> result;
    for (char c : test)
        result.push_back(std::string(1, c));
    return result;
}

const std::vector<std::string> Util::getToyExample1_Ans1() 
{
    std::string ans = "HHHLLLLLL";
    std::vector<std::string> result;
    for (char c : ans)
        result.push_back(std::string(1, c));
    return result;
}

const std::vector<std::string> Util::getToyExample1_Test2()
{
    std::string test = "GGCA";
    std::vector<std::string> result;
    for (char c : test)
        result.push_back(std::string(1, c));
    return result;
}

const std::vector<std::string> Util::getToyExample1_Ans2()
{
    std::string ans = "HHHL";
    std::vector<std::string> result;
    for (char c : ans)
        result.push_back(std::string(1, c));
    return result;  
}