// data generator
#include "util.h"
#include <iostream>
#include <string>
#include<stdlib.h>
using namespace std;

// num1 -- number of states
// num2 -- number of observations inside state
void generate(num1, num2){
    // generate states
    Table trans_prob;
    for(int i = 0; i<num1; i++)
    {
        String name_from = "State" + i;
        float probs[num1];
        float sum = 0.0;
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
        trans_prob.insert({name_from, {}});
        for (int j = 0; j<num1; j++)
        {
            String name_to = "State" + j;
            trans_prob[name_from].insert({name_to,probs[j]});
        }
    }
    cout << trans_prob;
}

int main(){
    generate(10,5);
    return 0;
}