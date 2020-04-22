// data generator
#include "util.h"
#include <iostream>
#include <string>
#include<stdlib.h>
using namespace std;

// num1 -- number of states
// num2 -- number of observations inside state
void generate(int num1, int num2){
    // generate states
    Table trans_prob;
    for(int i = 0; i<num1; i++)
    {
        string name_from = "State" + i;
        double probs[num1];
        double sum = 0.0;
        // Create array of random numbers
        for(int i=1;i<=num1;i++)
        {
            probs[i]=(rand() % 100) + 1;
            sum += probs[i];
        }
        // Divide each element to sum 
        for(int i=1;i<=num1;i++)
        {
            probs[i] /= sum;
        }
        map<string,double> states;
        trans_prob.insert(make_pair(name_from, states));
        for (int j = 0; j<num1; j++)
        {
            string name_to = "State" + j;
            
            trans_prob[name_from].insert(make_pair(name_to,probs[j]));
        }
    }
    for (auto& t : trans_prob)
    cout << t.first << " " 
         << t.second.first << " " 
         << t.second.second << "\n";
}

int main(){
    generate(10,5);
    return 0;
}