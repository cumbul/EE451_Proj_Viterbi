#include "util.h"
#include <iostream>
#include <string>
#include <json/value.h>
#include <fstream>

/****************************
Read data from JSON 
****************************/
// Read JSON data we generated from generator.py and add it into map
/*HMM Util::*/ void getData()
{
	//Table test_trans_prob;
	map<string, double> init_prob;
    std::ifstream init_file("data_init.json", std::ifstream::binary);
    init_file >> init_prob;
    cout<<init_prob;
    // Table test_obs_prob{
    // map<string, double> init_prob{
    // return HMM(test_trans_prob, test_obs_prob, init_prob);
}

int main()
{
	getData();
}