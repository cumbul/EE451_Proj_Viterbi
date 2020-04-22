// data generator
#include "util.h"
#include <iostream>
#include <string>
#include<stdlib.h>
#include <cstdlib>
using namespace std;

// num1 -- number of states
// num2 -- number of observations inside state
void generate(int num1, int num2){
    // generate states
    Table trans_prob;
    for(int i = 0; i<num1; i++)
    {

        string name_from = "State" + to_string(i);
        double probs[num1];
        double sum = 0.0;
        // Create array of random numbers
        for(int j=0;j<num1;j++)
        {
            srand(((i+10)*(j+1)) * 2654435789 + i);
            probs[j]=( rand() % 100) + 1;
            sum += probs[j];
        }
        // Divide each element to sum 
        for(int i=0;i<num1;i++)
        {
            probs[i] /= sum;
        }
        map<string,double> states;
        trans_prob.insert(make_pair(name_from, states));
        for (int j = 0; j<num1; j++)
        {
            string name_to = "State" + to_string(j);
            
            trans_prob[name_from].insert(make_pair(name_to,probs[j]));
        }
    }

    // generate observations
    Table obs_prob;
    for(int i = 0; i<num1; i++)
    {

        string name_from = "State" + to_string(i);
        double probs[num2];
        double sum = 0.0;
        // Create array of random numbers
        for(int j=0;j<num2;j++)
        {
            srand(((i+10)*(j+1)) * 90234859028 + i);
            probs[j]=( rand() % 100) + 1;
            sum += probs[j];
        }
        // Divide each element to sum 
        for(int i=0;i<num2;i++)
        {
            probs[i] /= sum;
        }
        map<string,double> states;
        obs_prob.insert(make_pair(name_from, states));
        for (int j = 0; j<num2; j++)
        {
            string name_to = "Obs" + to_string(j);
            
            obs_prob[name_from].insert(make_pair(name_to,probs[j]));
        }
    }
    // generate initial value
    map<string, double> init_prob;
    double probs[num1];
    double sum = 0.0;
    // Create array of random numbers
    for(int j=0;j<num1;j++)
    {
        srand((10*(j+1)) * 90234859028 + j);
        probs[j]=( rand() % 100) + 1;
        sum += probs[j];
    }
    // Divide each element to sum 
    for(int i=0;i<num1;i++)
    {
        probs[i] /= sum;
    }

    for(int i = 0; i<num1; i++)
    {
        string name_from = "State" + to_string(i);
        init_prob.insert(make_pair(name_from, probs[i]));
    }
    // ************
    // printing 
    // ************

    for (map<string, map<string, double>>::iterator it = trans_prob.begin(); it != trans_prob.end(); ++it){
        cout << it->first << " " << '\n';
        for (map<string, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            cout << it2->first << " " << it2->second << '\n';
    }

    for (map<string, map<string, double>>::iterator it = obs_prob.begin(); it != obs_prob.end(); ++it){
        cout << it->first << " " << '\n';
        for (map<string, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            cout << it2->first << " " << it2->second << '\n';
    }

    for (map<string, double>::iterator it = init_prob.begin(); it != init_prob.end(); ++it){
        cout << it->first << " " << it->second << '\n';
    }        
}

int main(){
    generate(10,5);
    return 0;
}