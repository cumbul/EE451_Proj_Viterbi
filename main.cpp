#include <iostream>

#include "viterbi.h"
#include "parallel_viterbi.h"
#include "util.h"
#include <map>

using namespace std;
//typedef map<string, map<string, float> > Table; 

int main()
{
    Viterbi viterbi(Util::getToyExample1());

    cout << "Vanilla Viterbi Test:" << endl;

    vector<string> input1 = Util::getToyExample1_Test1();
    cout << "Input 1: ";
    for (string s: input1) cout << s << " ";
    cout << endl;
    vector<string> my_ans = viterbi.solve(input1);
    cout << "My answer: ";
    for (string s: my_ans) cout << s << " ";
    cout << endl;
    vector<string> true_ans = Util::getToyExample1_Ans1();
    cout << "True answer: ";
    for (string s: true_ans) cout << s << " ";
    cout << endl;

    vector<string> input2 = Util::getToyExample1_Test2();
    cout << "Input 2: ";
    for (string s: input2) cout << s << " ";
    cout << endl;
    my_ans = viterbi.solve(input2);
    cout << "My answer: ";
    for (string s: my_ans) cout << s << " ";
    cout << endl; 
    true_ans = Util::getToyExample1_Ans2();
    cout << "True answer: ";
    for (string s: true_ans) cout << s << " ";
    cout << endl;

    viterbi = Viterbi(Util::getToyExample2());
    vector<string> input3 = Util::getToyExample2_Test();
    cout << "Input 3: ";
    for (string s: input3) cout << s << " ";
    cout << endl;
    my_ans = viterbi.solve(input3);
    cout << "My answer: ";
    for (string s: my_ans) cout << s << " ";
    cout << endl; 
    true_ans = Util::getToyExample2_Ans();
    cout << "True answer: ";
    for (string s: true_ans) cout << s << " ";
    cout << endl;

    LTDPViterbi ltdp = LTDPViterbi(Util::getToyExample2(), 2);
    my_ans = ltdp.solve(Util::getToyExample2_Test());
    cout << "LTDP Viterbi: ";
    for (string s: my_ans) cout << s << " ";
    cout << endl;    
}