//
// Created by kamilla on 03.04.19.
//

#ifndef FRESS_CPP_FRESS_H
#define FRESS_CPP_FRESS_H


#include <list>
#include <valarray>
#include <algorithm>
#include <random>
#include <tuple>
#include <vector>
#include <map>
#include <math.h>


class Protein{
public:
    std::vector <int> sequence;
    int number_of_iterations;
    int E;
    int min_E;
    float T;
    std::vector <double>  probabilities={};
    //std::pair <int, std::tuple<int, int>> node;
    //std::pair <int, int> coordinate;
    std::vector <std::pair <int, int>> conformation;
    std::vector<int> conformation_int;
    std:: vector <std:: list <std:: pair <int, int>>> results;
    std:: map <int, std::vector < std::pair <int, int> >> map_of_contacts;
    std:: map <int, std::pair<int, int>> map_int_to_coordinate;
    std:: map <std::pair<int, int>, int> map_coordinate_to_int;



    Protein();
    Protein(std::vector <int> sequence_input, int number_of_iterations=5000000 );

    void find_minimum();
    void calculate_probabilities_for_l(int lmin = 2, int lmax = 12);
    int count_contacts();
    void regrowth_middle(int l, int start_position);
    int distance( std:: pair <int, int> point1, std:: pair <int, int>point2   );



};

int count_contacts_breaked(std::vector <int> &sequence, std::vector <std::pair <int, int>> &conformation, std:: map <int, std::vector < std::pair <int, int> >> &map_of_contacts, std:: map <std::pair<int, int>, int> &map_coordinate_to_int     );







#endif //FRESS_CPP_FRESS_H
