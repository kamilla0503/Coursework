//
// Created by kamilla on 03.04.19.
//

#include "fress.h"
#include <iostream>

Protein::Protein(){};

Protein::Protein(std::valarray<int> sequence_input, int number_of_iterations_input) {
    sequence=sequence_input;
    number_of_iterations=number_of_iterations_input;
    calculate_probabilities_for_l();
    int length = sequence.size();
    conformation.push_back(std::make_pair(0, 0));
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

int Protein::count_contacts(){
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
    /**for ( std::pair <int, int> c : not_topological    ){
        std::cout << c.first  << " " << c.second << " ";
    }
    std::cout<< std:: endl;**/

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

        /**
       std::cout << "size of not topological : " << not_topological.size() << std::endl;
       for ( std::pair <int, int> q : not_topological    ){
           std::cout << q.first  << " " << q.second << " ";
       }
       std::cout<< std:: endl; **/


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

}

int Protein::distance( std:: pair <int, int> point1, std:: pair <int, int>point2   ){


    return abs(point2.first-point1.first) + abs(point2.second-point1.second);

}


void Protein::regrowth_middle(int l, int start_position){

    static std::list <std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
    std::pair <int, int> new_point ;
    int end_position = start_position+l-1;

    std::list <std::pair <int, int>> first_moves ;

    std::list <std::pair <int, int>> removed_part, copy, temporary_part;
    std::list<std::pair <int, int>>::iterator range_begin = conformation.begin();
    std::list<std::pair <int, int>>::iterator range_end = conformation.begin();
    std::advance(range_begin, start_position);
    std::advance(range_end, end_position+1);
    removed_part.splice(removed_part.begin(), conformation, range_begin, range_end   );
    copy=removed_part; //чтобы в случае неудачи вернуть предыдущую конформацию
    temporary_part = conformation;
   // std::list<std::pair <int, int>>::iterator range_begin_c = conformation.begin();
   // std::list<std::pair <int, int>>::iterator range_end_c = conformation.begin();
    range_begin=temporary_part.begin();
    //range_end=temporary_part.begin();
    std::advance(range_begin, end_position+1); //вроде на следующем после удаленного куска
    //std::advance(range_end, end_position+1);
    temporary_part.erase(range_begin, temporary_part.end()); // начальный кусок цепочки с выкинутой частью
    std::pair <int, int> finish_point = temporary_part.back();
    range_begin=temporary_part.begin();





    //std::cout << " l " << l << " size " << removed_part.size() << std::endl;




    for ( std::pair <int, int> step : steps ){
        new_point=std::make_pair(temporary_part.back().first+step.first, temporary_part.back().second+step.second);
        if(new_point==removed_part.front() ){
            continue;
        }
        if (std::find(conformation.begin(), conformation.end(), new_point)==conformation.end() &&(new_point, current_conf[end_position+1])  )




    }








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


        regrowth_middle(l, start_position);





    }
    std::cout << std::endl;

};