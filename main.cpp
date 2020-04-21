#include <iostream>

#include "viterbi.h"
#include "parallel_viterbi.h"
#include "util.h"
#include <map>

//typedef std::map<std::string, std::map<std::string, float> > Table; 

int main()
{
    Viterbi viterbi(Util::getToyExample1());

    std::cout << "Vanilla Viterbi Test:" << std::endl;

    std::vector<std::string> input1 = Util::getToyExample1_Test1();
    std::cout << "Input 1: ";
    for (std::string s: input1) std::cout << s << " ";
    std::cout << std::endl;
    std::vector<std::string> my_ans = viterbi.solve(input1);
    std::cout << "My answer: ";
    for (std::string s: my_ans) std::cout << s << " ";
    std::cout << std::endl;
    std::vector<std::string> true_ans = Util::getToyExample1_Ans1();
    std::cout << "True answer: ";
    for (std::string s: true_ans) std::cout << s << " ";
    std::cout << std::endl;

    std::vector<std::string> input2 = Util::getToyExample1_Test2();
    std::cout << "Input 2: ";
    for (std::string s: input2) std::cout << s << " ";
    std::cout << std::endl;
    my_ans = viterbi.solve(input2);
    std::cout << "My answer: ";
    for (std::string s: my_ans) std::cout << s << " ";
    std::cout << std::endl; 
    true_ans = Util::getToyExample1_Ans2();
    std::cout << "True answer: ";
    for (std::string s: true_ans) std::cout << s << " ";
    std::cout << std::endl;

    viterbi = Viterbi(Util::getToyExample2());
    std::vector<std::string> input3 = Util::getToyExample2_Test();
    std::cout << "Input 3: ";
    for (std::string s: input3) std::cout << s << " ";
    std::cout << std::endl;
    my_ans = viterbi.solve(input3);
    std::cout << "My answer: ";
    for (std::string s: my_ans) std::cout << s << " ";
    std::cout << std::endl; 
    true_ans = Util::getToyExample2_Ans();
    std::cout << "True answer: ";
    for (std::string s: true_ans) std::cout << s << " ";
    std::cout << std::endl;
}