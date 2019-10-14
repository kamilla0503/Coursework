//
// Created by kamilla on 13.10.2019.
//
#include "mcmc.h"
#include <iostream>
#include <fstream>


Lattice::Lattice() {};

Lattice::Lattice(int seq_size ) {
    lattice_side = seq_size+3;
    int l_s = seq_size+2; //чтобы было точно достатчно места при граничных условия для тора
    //создается словарь, в котором ключ - номер точки на квардатной решетке
    // значения - номера четырех соседних точек
    int x, y;
    div_t n;
    for (int i =0; i<(l_s+1)*(l_s+1); i++){
        //map_int_to_coordinate[i]=std::make_pair(i%(l_s+1), i /(l_s+1)  );
        //map_coordinate_to_int[map_int_to_coordinate[i]] = i;
        map_of_contacts_int[i] = {};
        map_of_contacts_int[i].push_back(i+1);
        map_of_contacts_int[i].push_back(i-1);
        map_of_contacts_int[i].push_back(i+l_s+1);
        map_of_contacts_int[i].push_back(i-l_s-1);
        n=div(i, l_s+1);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[i][1] = i+l_s;
            }
            if(x==l_s){
                map_of_contacts_int[i][0] = i-l_s;
            }
            if(y==0){
                map_of_contacts_int[i][3] = l_s*(l_s+1)+x;
            }
            if(y==l_s){
                map_of_contacts_int[i][2] = x;
            }

        }

    }

};

void Lattice::create_lattice(int seq_size = 900){
    lattice_side = seq_size+3;
    int l_s = seq_size+2; //чтобы было точно достатчно места при граничных условия для тора
    //создается словарь, в котором ключ - номер точки на квардатной решетке
    // значения - номера четырех соседних точек
    int x, y;
    div_t n;
    for (int i =0; i<(l_s+1)*(l_s+1); i++){
        //map_int_to_coordinate[i]=std::make_pair(i%(l_s+1), i /(l_s+1)  );
        //map_coordinate_to_int[map_int_to_coordinate[i]] = i;
        map_of_contacts_int[i] = {};
        map_of_contacts_int[i].push_back(i+1);
        map_of_contacts_int[i].push_back(i-1);
        map_of_contacts_int[i].push_back(i+l_s+1);
        map_of_contacts_int[i].push_back(i-l_s-1);
        n=div(i, l_s+1);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[i][1] = i+l_s;
            }
            if(x==l_s){
                map_of_contacts_int[i][0] = i-l_s;
            }
            if(y==0){
                map_of_contacts_int[i][3] = l_s*(l_s+1)+x;
            }
            if(y==l_s){
                map_of_contacts_int[i][2] = x;
            }

        }

    }

}

Protein::Protein(Sequence sequence_temp) {

    iterator=0;
    T=3.5;
    //change_T=false;
    lattice.create_lattice();
    sequence_template = sequence_temp;
    sequence = sequence_temp;
    for (int i=0; i<sequence_temp.size(); i++){
        conformation.push_back(i);
    }
    E = count_contacts();
    //min_E = E; //сохраняем минимальное найденное значение
    results_mean_E={};
}

Protein::Protein(){};



int Protein::count_contacts(){
    //считается число топологических контактов НН
    int hh = 0;
    int position;
    for (int i =1; i<sequence.size()-1; i++){
        for ( coord_t step : lattice.get_contacts(conformation[i]) ){
            if ( step!=conformation[i-1] && step!=conformation[i+1] && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){
                position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
                hh=hh+sequence[i]*sequence[position];
            }
        }
    }
    //концы конформации обрабатываются отдельно
    for ( coord_t step : lattice.get_contacts(conformation[0]) ) {
        if (step != conformation[1] &&
            std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {
            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
            hh = hh + sequence.front() * sequence[position];
        }
    }
    for ( coord_t step : lattice.get_contacts(conformation.back()) ){
        if ( step!=conformation[conformation.size()-2]  && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){
            position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
            hh=hh+sequence.back()*sequence[position];
        }
    }
    return  (-1*div(hh, 2).quot);
}




int Protein::dissected(Sequence &sequence1, Conformation &conformation1) {
    int hh = 0;
    int position;
    for (int i = 1; i < sequence1.size() - 1; i++) {
        for (coord_t step : lattice.get_contacts(conformation1[i])) {
            if (step != conformation1[i - 1] && step != conformation1[i + 1] &&
                std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
                position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
                hh = hh + sequence1[i] * sequence1[position];
            }
        }
    }
    for (coord_t step : lattice.get_contacts(conformation1[0])) {
        if (step != conformation1[1] &&
            std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
            position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
            hh = hh + sequence1[0] * sequence1[position];
        }
    }
    for (coord_t step : lattice.get_contacts(conformation1.back())) {
        if (step != conformation1[conformation1.size() - 2] &&
            std::find(conformation1.begin(), conformation1.end(), step) != conformation1.end()) {
            position = std::distance(conformation1.begin(), find(conformation1.begin(), conformation1.end(), step));
            hh = hh + sequence1.back() * sequence1[position];
        }
    }
    return (-1 * div(hh, 2).quot);


}

float Protein::MC_for_E(float temperature, float fugacity, int num_steps){
    float accumulate_E = 0.0;
    float accumulate_E_sq = 0.0;
    int r = 0;
    int a = 0;
    std::vector<int> lengths = {};
    std::vector<float > Es = {};
    int i=0;
    float q, w;
    Conformation new_conformation;
    Sequence new_sequence;
    int new_E = 0;
    float p_metropolis = 0.0;
    int rand_path;
    coord_t new_point;

    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);


    std::default_random_engine generators1(std::random_device{}() );
    std::uniform_int_distribution<int> distribution1(0, 3); //temp
    std::ofstream forl;          // поток для записи
    forl.open("leng_over_t_"+std::to_string(temperature)+".txt" );
    std::cout << "leng_over_t_"+std::to_string(temperature)+".txt" << std::endl;
    while(i<num_steps){
        w = distribution(generator);

        if(w<0.5){
            if(conformation.size()>1){
                new_conformation = conformation;
                new_sequence = sequence;
                new_conformation.pop_back();
                new_sequence.pop_back();
                new_E = dissected(new_sequence, new_conformation);
                p_metropolis = std::min(1.0, 1.0/exp(-(-new_E+E-fugacity)/temperature)/4.0 );
                q = distribution(generator);
                if(q<p_metropolis){
                    conformation = new_conformation;
                    sequence = new_sequence;
                    E = new_E;
                    r = r+1;
                }
            }
        }


        else {
            rand_path = distribution1(generators1);

            new_point = lattice.get_contacts(conformation.back())[rand_path];
            //std:: cout << rand_path << " " << new_point << std::endl;
            if(std::find(conformation.begin(), conformation.end(), new_point)==conformation.end()){
                new_conformation = conformation;
                new_sequence = sequence;
                new_conformation.push_back(new_point);
                new_sequence.push_back(sequence_template[sequence.size() % sequence_template.size()]);
                new_E = dissected(new_sequence, new_conformation);
                p_metropolis = std::min(1.0, exp(-(-new_E+E-fugacity)/temperature)*4.0 );
                q = distribution(generator);
                if(q<p_metropolis){
                    conformation = new_conformation;
                    sequence = new_sequence;
                    E = new_E;
                    a = a+1;
                }

            }

        }



        accumulate_E = accumulate_E + 1.0*(E-fugacity*sequence.size())/num_steps;

        forl<< sequence.size() << " ";


        i=i+1;
    }


    forl<<std::endl;
    forl.close();
    //std::cout << accumulate_E<< std:: endl;
    return accumulate_E;
    
}