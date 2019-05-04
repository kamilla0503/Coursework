//
// Created by kamilla on 03.04.19.
//

#include "fress.h"
#include <iostream>
#include <fstream>
Protein::Protein(){};

Protein::Protein(std::vector<int> sequence_input, int number_of_iterations_input) {
    sequence=sequence_input;
    number_of_iterations=number_of_iterations_input;
    T=3.5;
    calculate_probabilities_for_l();
    int length = sequence.size();

    int lattice_size = (length+2)*(length+2);
    int l_s = length+2;

  /**  for (int i=1; i<l_s ; i++) {
        for (int j=1; j<l_s+1 ; j++){
            map_of_contcts[i*l_s +j] = { std::make_pair( j-1, i+1),std::make_pair(j-1, i-1 ), std::make_pair(j, i), std::make_pair( j-2, i)         };

        }

        //map_of_contcts[]


    }

    **/

    /**
    for (int j=1; j<l_s ; j++){
      //mb ok
        map_of_contcts[ j] = { std::make_pair(j, 1),std::make_pair(j, l_s), std::make_pair(j+1, 0), std::make_pair(j-1,0)         };
//!!!!
        map_of_contcts[l_s*(l_s-1)+ j] =  { std::make_pair(l_s-1+ j, 0),std::make_pair(l_s-1+ j, l_s-1), std::make_pair(l_s*(l_s-1)+ j+1, l_s), std::make_pair(l_s*(l_s-1)+ j-1,l_s)         };

        map_of_contcts[l_s*j+1] = { std::make_pair( 0, j+1),std::make_pair( 0, j-1), std::make_pair(1, j), std::make_pair(l_s, j)         };

        //map_of_contcts[l_s*j] = { std::make_pair( 0, j+1),std::make_pair( 0, j-1), std::make_pair(1, j), std::make_pair(l_s, j)         };

        map_of_contcts[l_s*(1+j) +1] = { std::make_pair( l_s, j+1),std::make_pair( l_s, j-1), std::make_pair(l_s-1, j), std::make_pair(0, j)         };

    }**/


    for (int i=0; i<=l_s ; i++) {
        for (int j=0; j<=l_s ; j++){
            map_of_contacts[i*(l_s+1) +j] = { std::make_pair( j-1, i),std::make_pair(j+1, i ), std::make_pair(j, i-1), std::make_pair( j, i+1)         };
            //for (std::pair<int, int> c : map_of_contacts[i*l_s+j] ){
            for (int c=0; c<  map_of_contacts[i*l_s +j].size(); c++ ){
                if(map_of_contacts[i*(l_s+1) +j][c].first>l_s){
                    map_of_contacts[i*(l_s+1) +j][c].first = 0;
                }
                if(map_of_contacts[i*(l_s+1) +j][c].second>l_s){
                    map_of_contacts[i*(l_s+1) +j][c].second = 0;
                }
                if(map_of_contacts[i*(l_s+1) +j][c].second<0){
                    map_of_contacts[i*(l_s+1) +j][c].second = l_s;
                }
                if(map_of_contacts[i*(l_s+1) +j][c].first<0){
                    map_of_contacts[i*(l_s+1) +j][c].first = l_s;
                }
            }

        }


    }

    for (int i =0; i<(l_s+1)*(l_s+1); i++){
        map_int_to_coordinate[i]=std::make_pair(i%(l_s+1), i /(l_s+1)  );
        map_coordinate_to_int[map_int_to_coordinate[i]] = i;
    }
/**
    for (auto c : map_int_to_coordinate) {

        std:: cout << c.first << " " << c.second.first << " " << c.second.second << std::endl;
    }
**/
    for (auto c : map_coordinate_to_int) {

        std:: cout << c.first.first << " " << c.first.second << " " << c.second  << std::endl;
    }


    std::ofstream out1;
    out1.open("contact_map");

    for (auto c : map_of_contacts){
        out1 << c.first << " ";
        for (auto d : c.second){
            out1 << d.first << " " << d.second << "   ";
        }

        out1 << std::endl;
    }



    conformation.push_back(std::make_pair(0, 0));
    conformation_int.push_back(0);
    std::pair <int, int> new_coordinate;
    for (int i=1; i<length; i++){
        if (i%8>0 && i%8<4) {
            new_coordinate=std::make_pair(conformation.back().first, conformation.back().second+1);
        }
        else if (i%8>4 && i%8<=7){
            new_coordinate=std::make_pair(conformation.back().first, conformation.back().second-1);
        }
        else{
            new_coordinate=std::make_pair(conformation.back().first+1, conformation.back().second);
        }
        conformation.push_back(new_coordinate);
        conformation_int.push_back(map_coordinate_to_int[new_coordinate]);
    }

    std::ofstream out;
    out.open("coordinates_one_step");

    for (auto c : conformation) {

        out << c.first << " " << c.second << "   ";


        out << std::endl;


    }






    E = count_contacts();
    min_E = E;
    results={};

    /**std::cout << "start conformation size : " << conformation.size() << std::endl;
    for ( std::pair <int, int> c : conformation     ){
        std::cout << c.first  << " " << c.second << std::endl;
    }**/

}


void Protein::calculate_probabilities_for_l(int lmin, int lmax) {
    double one_over_l = 0.0;
    for (int i=lmin; i <=lmax; i++ ){
        one_over_l=one_over_l+1.0/i;
    }

   // std::cout << lmin << " " << lmax << std::endl;
    double  k = 1/one_over_l;
   // std::cout << "one over l " << one_over_l << std::endl;
    for (int i=lmin; i <=lmax; i++ ){
        probabilities.push_back(k/i);
    }

    for (int i =1; i<probabilities.size(); i++){
        probabilities[i]=probabilities[i]+probabilities[i-1];
    }


}
//something old
/**  int Protein::count_contacts(){
    int hh=0;
    static std::list <std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
    std:: list <std::pair <int, int>> not_topological;
    not_topological=conformation;
    int i =1;
    int position;
    std::list<std::pair <int, int>>::iterator range_begin = not_topological.begin();
    std::list<std::pair <int, int>>::iterator range_end = not_topological.begin();
    std::advance(range_begin,3);
    std::advance(range_end, conformation.size());
    not_topological.erase(range_begin, range_end);
    std::list<std::pair <int, int>>::iterator range_middle = not_topological.begin();
    std::advance(range_middle,1);
    not_topological.erase(range_middle);
    std::pair <int, int> new_point, new_point1;
    //std::cout << "size of not topological : " << not_topological.size() << std::endl;

   for (std::pair <int, int> c : conformation ) {

       if(c!=conformation.front() && c!=conformation.back()) {
           for ( std::pair <int, int> step : steps ) {
               new_point=std::make_pair( c.first+step.first, c.second+step.second );
               if ( std::find(conformation.begin(), conformation.end(), new_point) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point)==not_topological.end()    ) {
                   position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), new_point));
                hh=hh+sequence[i]*sequence[position];
               }
           }

           i=i+1;
           not_topological.clear();
           if(i!=sequence.size()-3 ) {
               not_topological = conformation;
               std::list<std::pair <int, int>>::iterator range_middle = not_topological.begin();
               std::list<std::pair <int, int>>::iterator range_end = not_topological.begin();
               std::advance(range_middle, i+2);
               std::advance(range_end, conformation.size());
               not_topological.erase(range_middle, range_end);
               new_point1=not_topological.back();
               not_topological.clear();
               not_topological.push_back(c);
               not_topological.push_back(new_point1);

           }
           else if(i==sequence.size()-3){
               new_point1=conformation.back();
               not_topological.push_back(c);
               not_topological.push_back(new_point1);
           }

       }



   }

    for ( std::pair <int, int> step : steps ){
        std::pair <int, int> new_point_begin = std::make_pair(conformation.front().first+step.first, conformation.front().second+step.second );
        std::pair <int, int> new_point_end= std::make_pair(conformation.back().first+step.first,conformation.back().second+step.second);
        not_topological.clear();
        not_topological=conformation;
        std::list<std::pair <int, int>>::iterator range_middle = not_topological.begin();
        std::list<std::pair <int, int>>::iterator range_end = not_topological.begin();
        std::advance(range_middle, 2);
        std::advance(range_end, conformation.size()-3);
        not_topological.erase(range_middle, range_end);
        not_topological.pop_front();
        not_topological.pop_back(); //теперь тут только второй и предпоследний элементы
        if (std::find(conformation.begin(), conformation.end(), new_point_begin) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point_begin)==not_topological.end()   ) {
            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), new_point_begin) );
            hh = hh + sequence[position]*sequence[0];
        }

        if (std::find(conformation.begin(), conformation.end(), new_point_end) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point_end)==not_topological.end()   ) {
            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), new_point_end) );
            hh = hh + sequence[position]*sequence[0];
        }

    }

    return (-1*div(hh, 2).quot);

}**/

/**
int count_contacts(std::vector<std::pair <int, int>> &conformation, std::valarray<int> &sequence){
    int hh = 0;
    int position;
    static std::valarray<std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
    std:: vector <std::pair <int, int>> not_topological = {};
    std::pair <int, int> new_point,new_point_begin, new_point_end;
    for (int i =1; i<sequence.size()-1; i++){
        not_topological.push_back(conformation [i-1]);
        not_topological.push_back(conformation [i+1]);
        for ( std::pair <int, int> step : steps ){
            new_point = std::make_pair( conformation[i].first+step.first, conformation [i].second+step.second );
            if ( std::find(conformation.begin(), conformation.end(), new_point) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point)==not_topological.end() ){
                position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), new_point));
                hh=hh+sequence[i]*sequence[position];
            }



        }

        not_topological={};


    }


    for ( std::pair <int, int> step : steps ) {
        new_point_begin = std::make_pair(conformation[0].first+step.first, conformation[0].second+step.second );
        new_point_end= std::make_pair(conformation.back().first+step.first,conformation.back().second+step.second);
        if(std::find(conformation.begin(), conformation.end(), new_point) !=conformation.end()  && new_point_begin!= conformation [1]) {
            position = std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), new_point_begin));

            hh=hh+sequence[position]*sequence[0];


        }


        if (std::find(conformation.begin(), conformation.end(), new_point_end) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point_end)==not_topological.end()   ) {
            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), new_point_end) );
            hh = hh + sequence[position]*sequence[sequence.size()-1];
        }






    }




    return  (-1*div(hh, 2).quot);

} **/



/**int Protein::count_contacts(){
    int hh = 0;
    int position;
    static std::valarray<std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
    std:: vector <std::pair <int, int>> not_topological = {};
    std::pair <int, int> new_point,new_point_begin, new_point_end;
    for (int i =1; i<sequence.size()-1; i++){
        not_topological.push_back(conformation [i-1]);
        not_topological.push_back(conformation [i+1]);
        for ( std::pair <int, int> step : steps ){
            new_point = std::make_pair( conformation[i].first+step.first, conformation [i].second+step.second );
            if ( std::find(conformation.begin(), conformation.end(), new_point) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point)==not_topological.end() ){
                position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), new_point));
                hh=hh+sequence[i]*sequence[position];
            }



        }

        not_topological={};


    }


    for ( std::pair <int, int> step : steps ) {
        new_point_begin = std::make_pair(conformation[0].first+step.first, conformation[0].second+step.second );
        new_point_end= std::make_pair(conformation.back().first+step.first,conformation.back().second+step.second);
        if(std::find(conformation.begin(), conformation.end(), new_point) !=conformation.end()  && new_point_begin!= conformation [1]) {
            position = std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), new_point_begin));

            hh=hh+sequence[position]*sequence[0];


        }


        if (std::find(conformation.begin(), conformation.end(), new_point_end) !=conformation.end() &&std::find(not_topological.begin(),not_topological.end(), new_point_end)==not_topological.end()   ) {
            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), new_point_end) );
            hh = hh + sequence[position]*sequence[sequence.size()-1];
        }






    }




    return  (-1*div(hh, 2).quot);

}
**/


int Protein::count_contacts(){
    int hh = 0;
    int position;
    //static std::valarray<std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
    //std:: vector <std::pair <int, int>> not_topological = {};
    std::pair <int, int> new_point,new_point_begin, new_point_end;
    for (int i =1; i<sequence.size()-1; i++){
        //not_topological.push_back(conformation [i-1]);
       // not_topological.push_back(conformation [i+1]);
        for ( std::pair <int, int> step : map_of_contacts[map_coordinate_to_int[conformation[i]]] ){
            //new_point = std::make_pair( conformation[i].first+step.first, conformation [i].second+step.second );
            if ( step!=conformation[i-1] && step!=conformation[i+1] && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){

                position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
                hh=hh+sequence[i]*sequence[position];
            }



        }

        //not_topological={};


    }

    for ( std::pair <int, int> step : map_of_contacts[map_coordinate_to_int[conformation[0]]] ) {
        if (step != conformation[1] &&
            std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {

            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
            hh = hh + sequence.front() * sequence[position];


        }
    }

    for ( std::pair <int, int> step : map_of_contacts[map_coordinate_to_int[conformation.back()]] ){
        if ( step!=conformation[conformation.size()-2]  && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){

            position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
            hh=hh+sequence.back()*sequence[position];


        }


    }


    return  (-1*div(hh, 2).quot);

}




/**
int Protein::distance( std:: pair <int, int> point1, std:: pair <int, int>point2   ){


    return abs(point2.first-point1.first) + abs(point2.second-point1.second);

}**/

int Protein::distance( std:: pair <int, int> point1, std:: pair <int, int>point2   ){
    int p1 = abs(point2.first-point1.first);
    int p2 = sequence.size()+3 -  abs(point2.first-point1.first);
    int part1 = std::min(p1, p2);
    p1 = abs(point2.second-point1.second);
    p2 = sequence.size()+3 -  abs(point2.second-point1.second);
    int part2 = std::min(p1, p2);

    return  part1 + part2;

}





void Protein::regrowth_middle(int l, int start_position){

    //static std::list <std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };

    int end_position = start_position+l-1;
    std::vector<std::pair <int, int>>  C_t, C_t_temp ;
    std::vector<int> seq_t, seq_t_temp;
    //std::copy(sequence.begin(), seq_t.begin()+ start_position, C_t_temp.begin());
    C_t_temp.resize(start_position);
    //C_t_temp.resize(sequence.size());
    std::copy(conformation.begin(), conformation.begin()+ start_position, C_t_temp.begin());
    C_t = C_t_temp;
    //C_t_temp.clear();

    C_t_temp.resize(sequence.size()-start_position-l);
    std::copy(conformation.begin()+end_position+1, conformation.end()+ start_position, C_t_temp.begin());

    C_t.reserve(C_t.size()+C_t_temp.size());
    C_t.insert(C_t.end(), C_t_temp.begin(), C_t_temp.end());
    C_t_temp.clear();
    //это удаление куска конформации с start_position до end_position


    seq_t_temp.resize(start_position);
    //seq_t_temp.resize(sequence.size());
    std::copy(sequence.begin(), sequence.begin()+ start_position, seq_t_temp.begin());
    seq_t = seq_t_temp;
    seq_t_temp.clear();
    std::copy(sequence.begin()+end_position+1, sequence.end()+ start_position, seq_t_temp.begin());
    seq_t.reserve(seq_t.size()+seq_t_temp.size() );
    seq_t.insert(seq_t.end(), seq_t_temp.begin(), seq_t_temp.end());
    seq_t_temp.clear();
    //это удаление куска последовательности

    int current_energy = count_contacts_breaked(seq_t, C_t, map_of_contacts, map_coordinate_to_int);
    int temp_e;
    static std::vector <std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0, 0));
   // std::vector<std::pair <float , int>> first_moves ;
    std::vector<int> energies;
    energies.resize(4, 0);
    std::pair<int, int> point;
    float sum_probabilities=0.0;
    //for ( std::pair <int, int> step : map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]] ){
    for (int i =0; i<4; i++ ){
       // sum_probabilities=0.0;
        std:: cout << " check map"<< map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i].first << " " << map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i].second<< std:: endl;
        if (map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i] == conformation[start_position]){
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        else  if( std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i]   )==C_t.end() &&distance(map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i], conformation[end_position+1]) <=abs(end_position+1-start_position)  ){


            // Лучше потом переделать функцию для энергии
            C_t.insert(C_t.begin()+start_position+1,map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i]  );

            temp_e = count_contacts_breaked(sequence, C_t, map_of_contacts,map_coordinate_to_int);
            energies[i] = temp_e;
            probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

            C_t.erase(C_t.begin()+start_position);
        }
        else{
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }

        sum_probabilities=sum_probabilities+probabilities_to_move[i].first;



    }


    if(sum_probabilities==0.0){
        std:: cout << " fail to rebuild" << std::endl;
        return;

    }


   // std::cout << "sum for notm " << sum_probabilities <<std::endl;

    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


    }



    std::cout << "first move " << std:: endl;
    for (int i=0; i<4; i++){

        std::cout << probabilities_to_move[i].first <<" "<< probabilities_to_move[i].second << " " << std::endl;



        //std::cout << probabilities_to_move[i].first <<" "<< probabilities_to_move[i].second << " " <<   map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i].first << " " <<  map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i].second<< std::endl;

    }




    for (int i =1; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;


    }


    std::cout << "first move " << std:: endl;
    for (int i=0; i<4; i++){

        std::cout << probabilities_to_move[i].first <<" "<< probabilities_to_move[i].second << " " << std::endl;



        //std::cout << probabilities_to_move[i].first <<" "<< probabilities_to_move[i].second << " " <<   map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i].first << " " <<  map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i].second<< std::endl;

    }



    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator);
    std::cout << " q = " << q << std:: endl;
    for (int i =0; i<4; i++ ){

        if (q<probabilities_to_move[i].first){
            C_t.insert(C_t.begin()+start_position+1,map_of_contacts[map_coordinate_to_int[ C_t[start_position-1]  ]][probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];

            seq_t.insert(seq_t.begin()+start_position+1, sequence[start_position]);

            std:: cout << " choosed " << i << " " << std:: endl;
            std:: cout << "check size " << C_t.size() << std:: endl;
            break;

        }

    }



   // sum_probabilities = 0.0;


    for (int t = start_position+1; t< end_position+1; t++){
        sum_probabilities = 0.0;


        for (int i =0; i<4; i++ ){

            if(std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t[t-1]]][i])!=C_t.end() ||  distance( map_of_contacts[map_coordinate_to_int[C_t[t-1]]][i] ,conformation[end_position+1] )>abs(end_position+1-t)){

                probabilities_to_move[i] = std::make_pair(0.0, i);



            }
            else {
                //точка подходит

                C_t.insert(C_t.begin()+t+1,map_of_contacts[map_coordinate_to_int[C_t[t-1]]][i]  );

                temp_e = count_contacts_breaked(sequence, C_t, map_of_contacts,map_coordinate_to_int);
                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

                C_t.erase(C_t.begin()+t);





            }





            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;

        }



        if(sum_probabilities==0.0){
            std:: cout << " fail to rebuild" << std::endl;
            return;

        }

        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


        }
        for (int i =1; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].second;


        }




        q = distribution(generator);

        for (int i =0; i<4; i++ ){

            if (q<probabilities_to_move[i].first){
                C_t.insert(C_t.begin()+t+1,map_of_contacts[map_coordinate_to_int[ C_t[t-1]  ]][probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];

                seq_t.insert(seq_t.begin()+t+1, sequence[t]);



            }

        }










    }



// проверка по сути не нужна, но мне так спокойнее

    if( C_t.size()==sequence.size()  ){

        E = current_energy;
        conformation=C_t;

        std::cout << " check energy " << current_energy << std::endl;
        std::ofstream out;
        out.open("coordinates_one_step");

        for (auto c : conformation) {

            out << c.first << " " << c.second << "   ";


            out << std::endl;


        }


    }
    else{
    std::cout << " fail int the end " << std:: endl;
    return ;
    }

    sum_probabilities = 0.0;









}







void Protein::find_minimum() {

    double q;
    int l, start_position;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for (int iteration =1; iteration<number_of_iterations+1; iteration++) {


        q = distribution(generator);
        //std::cout<< q << " " ;

        for (int i = 0; i < probabilities.size(); i++) {
            if (q < probabilities[i]) {

                l =i+2;
                break;
            }


        }


        std::default_random_engine generators1;
        std::uniform_int_distribution<int> distribution1(0,sequence.size()-l-1);
        start_position = distribution1(generators1);


        //std:: cout << l << " " ;





    }
    std::cout << std::endl;

};

int count_contacts_breaked(std::vector <int> &sequence, std::vector <std::pair <int, int>> &conformation, std:: map <int, std::vector < std::pair <int, int> >> &map_of_contacts, std:: map <std::pair<int, int>, int> &map_coordinate_to_int    ) {
    int hh = 0;
    int position;
//static std::valarray<std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
//std:: vector <std::pair <int, int>> not_topological = {};
    std::pair<int, int> new_point, new_point_begin, new_point_end;
    for (int i = 1; i < sequence.size() - 1; i++) {
//not_topological.push_back(conformation [i-1]);
// not_topological.push_back(conformation [i+1]);
        for (std::pair<int, int> step : map_of_contacts[map_coordinate_to_int[conformation[i]]]) {
//new_point = std::make_pair( conformation[i].first+step.first, conformation [i].second+step.second );
            if (step != conformation[i - 1] && step != conformation[i + 1] &&
                std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {

                position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
                hh = hh + sequence[i] * sequence[position];
            }


        }

//not_topological={};


    }

    for (std::pair<int, int> step : map_of_contacts[map_coordinate_to_int[conformation[0]]]) {
        if (step != conformation[1] &&
            std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {

            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
            hh = hh + sequence[0] * sequence[position];


        }
    }

    for (std::pair<int, int> step : map_of_contacts[map_coordinate_to_int[conformation.back()]]) {
        if (step != conformation[conformation.size() - 2] &&
            std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {

            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
            hh = hh + sequence.back() * sequence[position];


        }






    }


    return (-1 * div(hh, 2).quot);

}