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


class Protein{
public:
    std::valarray <int> sequence;
    int number_of_iterations;
    int E;
    int min_E;
    std::vector <double>  probabilities={};
    //std::pair <int, std::tuple<int, int>> node;
    //std::pair <int, int> coordinate;
    std::vector <std::pair <int, int>> conformation;
    std:: vector <std:: list <std:: pair <int, int>>> results;
    Protein();
    Protein(std::valarray <int> sequence_input, int number_of_iterations=5000000 );

    void find_minimum();
    void calculate_probabilities_for_l(int lmin = 2, int lmax = 12);
    int count_contacts();
    void regrowth_middle(int l, int start_position);




};







#endif //FRESS_CPP_FRESS_H
