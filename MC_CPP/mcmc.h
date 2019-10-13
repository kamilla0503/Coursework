//
// Created by kamilla on 13.10.2019.
//
#ifndef MC_CPP_MCMC_H
#define MC_CPP_MCMC_H
#include <list>
#include <valarray>
#include <algorithm>
#include <random>
#include <tuple>
#include <vector>
#include <map>
#include <math.h>
typedef int coord_t;
class Lattice{
public:
    int ndim() { return 2; }
    int ndim2() {return 4;}
    int lattice_side;
    Lattice();
    Lattice(int seq_size );
    void create_lattice(int seq_size);
    //int distance_lattice(coord_t point1, coord_t point2);
    inline std::vector<coord_t> get_contacts(coord_t point){
        return map_of_contacts_int[point];
    }
private:
    std::map<int, std::vector<int>> map_of_contacts_int;
};
class Protein{
public:
    typedef std::vector<int> Sequence_t, Sequence;
    typedef std::vector<coord_t> Conformation_t, Conformation;
    Sequence sequence;
    Sequence sequence_template;
    int E; // энергия текущей конформации
    //int min_E;
    float T; //температура
    //bool change_T; //чтобы корректно менять температуру каждый 50 000 шагов
    long int iterator; //для подсчета числа шагов

    Conformation conformation;
    //std::vector<Conformation> results;
    std::vector<float> results_mean_E;
    Lattice lattice;
    Protein();
    Protein(Sequence sequence_temp);
    int count_contacts();
    int dissected(Sequence_t &sequence, Conformation_t &conformation);
    float MC_for_E(float temperature, float fugacity, int num_steps);

private:

};






#endif //MC_CPP_MCMC_H
