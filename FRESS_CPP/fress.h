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
   // int number_of_iterations;
    int E; // энергия текущей конформации
    int min_E;
    float T; //температура
    bool change_T; //чтобы корректно менять температуру каждый 50 000 шагов
    long int iterator; //для подсчета числа шагов
    std::vector <double>  probabilities={}; //вектор вероятностей для выбора длины, вероятность пропорциональна 1/l
    std::vector <std::pair <int, int>> conformation;
    std::vector<int> conformation_int;
    std:: vector <std:: vector <std:: pair <int, int>>> results;
    std:: map <int, std::vector < std::pair <int, int> >> map_of_contacts;
    std:: map <int, std::pair<int, int>> map_int_to_coordinate;
    std:: map <std::pair<int, int>, int> map_coordinate_to_int; //словарь для хранения всех соседей



    Protein();
    Protein(std::vector <int> sequence_input  );

    void find_minimum();
    void calculate_probabilities_for_l(int lmin = 2, int lmax = 12);
    int count_contacts();
    void regrowth_middle(int l, int start_position);
    void regrowth_start(int l );
    void regrowth_end(int l );
    int distance( std:: pair <int, int> point1, std:: pair <int, int>point2   );



};

int count_contacts_breaked(std::vector <int> &sequence, std::vector <std::pair <int, int>> &conformation, std:: map <int, std::vector < std::pair <int, int> >> &map_of_contacts, std:: map <std::pair<int, int>, int> &map_coordinate_to_int     );







#endif //FRESS_CPP_FRESS_H
