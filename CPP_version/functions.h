//
// Created by kamilla on 09.12.18.
//

#ifndef CPP_VERSION_FUNCTIONS_H
#define CPP_VERSION_FUNCTIONS_H

#include<vector>
#include<map>
#include <tuple>
#include <iostream>
#include <algorithm>
using namespace std;

vector <int> vector_for_distance(vector<tuple<int, int>> saw );
int distance_between_saws(vector<tuple<int, int>> saw1, vector<tuple<int, int>> saw2);
vector<vector<tuple<int, int>>>  filter_conformations(vector<vector<tuple<int, int>>>  saws );
vector<vector<tuple<int, int>>> get_all_conformations(int length);



#endif //CPP_VERSION_FUNCTIONS_H
