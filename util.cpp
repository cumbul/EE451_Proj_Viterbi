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

HMM::HMM(const Table& transition_prob, const Table& observation_prob, const map<string, double>& init_prob)
{
    this->trans_prob = transition_prob;
    this->obs_prob = observation_prob;
    this->init_prob = init_prob;
}

/*
Get transition probability from previous state to the next state.
*/
double HMM::get_trans_prob(const string& previous, const string& next)
{
    //If the input not found in table, show warning. To be deleted after debugging.
    if (this->trans_prob.find(previous) == this->trans_prob.end())
        cout << "Warning: previous argument '" << previous << "' in HMM::get_trans_prob() does not exist. Return 0 by default." << endl; 
    else if (this->trans_prob[previous].find(next) == this->trans_prob[previous].end())
        cout << "Warning: next argument '" << next << "' in HMM::get_trans_prob() does not exist. Return 0 by default." << endl; 

    return this->trans_prob[previous][next];
}

/*
Get observation likelihood of the observation given the current state.
*/
double HMM::get_obs_prob(const string& state, const string& observation)
{
    //If the input not found in table, show warning. To be deleted after debugging.
    if (this->obs_prob.find(state) == this->obs_prob.end())
        cout << "Warning: state argument '" << state << "' in HMM::get_obs_prob() does not exist. Return 0 by default." << endl; 
    else if (this->obs_prob[state].find(observation) == this->obs_prob[state].end())
        cout << "Warning: observation argument '" << observation << "' in HMM::get_obs_prob() does not exist. Return 0 by default." << endl; 

    return this->obs_prob[state][observation];
}

double HMM::get_init_prob(const string& state)
{
    //If the input not found in table, show warning. To be deleted after debugging.
    if (this->init_prob.find(state) == this->init_prob.end())
        cout << "Warning: state argument '" << state << "' in HMM::get_init_prob() does not exist. Return 0 by default." << endl;
    return this->init_prob[state];
}

/*
Get a list of all possible states in HMM.
*/
vector<string> HMM::get_state_list()
{
    vector<string> result;
    for(Table::iterator it = this->trans_prob.begin(); it != this->trans_prob.end(); it++)
    {
        result.push_back(it->first);
    }
    return result;
}

/*
Get a list of all possible observation in HMM.
*/
vector<string> HMM::get_observation_list()
{
    vector<string> result;
    map<string, double> first_table = this->obs_prob.begin()->second;
    for(map<string, double>::iterator it = first_table.begin(); it != first_table.end(); it++)
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
    map<string, double> init_prob{
        {"H", 0.5},
        {"L", 0.5}
    };
    return HMM(test_trans_prob, test_obs_prob, init_prob);
}

const vector<string> Util::getToyExample1_Test1()
{
    string test = "GGCACTGAA";
    vector<string> result;
    for (char c : test)
        result.push_back(string(1, c));
    return result;
}

const vector<string> Util::getToyExample1_Ans1() 
{
    string ans = "HHHLLLLLL";
    vector<string> result;
    for (char c : ans)
        result.push_back(string(1, c));
    return result;
}

const vector<string> Util::getToyExample1_Test2()
{
    string test = "GGCA";
    vector<string> result;
    for (char c : test)
        result.push_back(string(1, c));
    return result;
}

const vector<string> Util::getToyExample1_Ans2()
{
    string ans = "HHHL";
    vector<string> result;
    for (char c : ans)
        result.push_back(string(1, c));
    return result;  
}

/****************************
Toy HMM 2 model for testing
****************************/
//return HMM model in page 155 of https://web.stanford.edu/~jurafsky/slp3/ed3book.pdf 
HMM Util::getToyExample2()
{
    Table test_trans_prob{
		{"NNP",
            {{"NNP", 0.3777}, {"MD", 0.0110}, {"VB", 0.0009}, {"JJ", 0.0084}, {"NN", 0.0584}, {"RB", 0.0090}, {"DT", 0.0025}}
        },
		{"MD",
            {{"NNP", 0.0008}, {"MD", 0.0002}, {"VB", 0.7968}, {"JJ", 0.0005}, {"NN", 0.0008}, {"RB", 0.1698}, {"DT", 0.0041}}
        },     
		{"VB",
            {{"NNP", 0.0322}, {"MD", 0.0005}, {"VB", 0.0050}, {"JJ", 0.0837}, {"NN", 0.0615}, {"RB", 0.0514}, {"DT", 0.2231}}
        },
		{"JJ",
            {{"NNP", 0.0366}, {"MD", 0.0004}, {"VB", 0.0001}, {"JJ", 0.0733}, {"NN", 0.4509}, {"RB", 0.0036}, {"DT", 0.0036}}
        },
		{"NN",
            {{"NNP", 0.0096}, {"MD", 0.0176}, {"VB", 0.0014}, {"JJ", 0.0086}, {"NN", 0.1216}, {"RB", 0.0177}, {"DT", 0.0068}}
        },
		{"RB",
            {{"NNP", 0.0068}, {"MD", 0.0102}, {"VB", 0.1011}, {"JJ", 0.1012}, {"NN", 0.0120}, {"RB", 0.0728}, {"DT", 0.0479}}
        },
		{"DT",
            {{"NNP", 0.1147}, {"MD", 0.0021}, {"VB", 0.0002}, {"JJ", 0.2157}, {"NN", 0.4744}, {"RB", 0.0102}, {"DT", 0.0017}}
        }
    };   
    Table test_obs_prob{
        {"NNP",
            {{"Janet", 0.000032}, {"will", 0}, {"back", 0}, {"the", 0.000048}, {"bill", 0}}
        },
        {"MD",
            {{"Janet", 0}, {"will", 0.308431}, {"back", 0}, {"the", 0}, {"bill", 0}}
        },
        {"VB",
            {{"Janet", 0}, {"will", 0.000028}, {"back", 0.000672}, {"the", 0}, {"bill", 0.000028}}
        },
        {"JJ",
            {{"Janet", 0}, {"will", 0}, {"back", 0.000340}, {"the", 0}, {"bill", 0}}
        },
        {"NN",
            {{"Janet", 0}, {"will", 0.000200}, {"back", 0.000223}, {"the", 0}, {"bill", 0.002337}}
        },
        {"RB",
            {{"Janet", 0}, {"will", 0}, {"back", 0.010446}, {"the", 0}, {"bill", 0}}
        },
        {"DT",
            {{"Janet", 0}, {"will", 0}, {"back", 0}, {"the", 0.506099}, {"bill", 0}}
        },
    };
    map<string, double> init_prob{
        {"NNP", 0.2767}, {"MD", 0.0006}, {"VB", 0.0031}, {"JJ", 0.0453}, {"NN", 0.0449}, {"RB", 0.0510}, {"DT", 0.2026}
    };
    return HMM(test_trans_prob, test_obs_prob, init_prob);
}

const vector<string> Util::getToyExample2_Test()
{
    vector<string> result;
    result.push_back("Janet");
    result.push_back("will");
    result.push_back("back");
    result.push_back("the");
    result.push_back("bill");
    return result; 
}

const vector<string> Util::getToyExample2_Ans()
{
    vector<string> result;
    result.push_back("NNP");
    result.push_back("MD");
    result.push_back("VB");
    result.push_back("DT");
    result.push_back("NN");
    return result;
}
