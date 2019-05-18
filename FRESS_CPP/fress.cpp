//
// Created by kamilla on 03.04.19.
//

#include "fress.h"
#include <iostream>
#include <fstream>
Protein::Protein(){};

Protein::Protein(std::vector<int> sequence_input ) {
    sequence=sequence_input;
    iterator=0;
    T=3.5;
    change_T=false;
    calculate_probabilities_for_l();
    int length = sequence.size();
    int lattice_size = (length+2)*(length+2);
    int l_s = length+2; //чтобы было точно достатчно места при граничных условия для тора


    //создается словарь, в котором ключ - номер точки на квардатной решетке
    // значения - координаты четырех соседних точек
    // создается пока неэффективным и некрасивым путем, но эта генерация происходит только один раз

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

// push_back or emplace_back?
    conformation.emplace_back(std::make_pair(0, 0));
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

    E = count_contacts();
    min_E = E; //сохраняем минимальное найденное значение
    results={};

}


void Protein::calculate_probabilities_for_l(int lmin, int lmax) {
    //создается вектор с вероятностями, с которыми будут выбираться l
    double one_over_l = 0.0;
    for (int i=lmin; i <=lmax; i++ ){
        one_over_l=one_over_l+1.0/i;
    }

    double  k = 1/one_over_l;
    for (int i=lmin; i <=lmax; i++ ){
        probabilities.push_back(k/i);
    }

    for (int i =1; i<probabilities.size(); i++){
        probabilities[i]=probabilities[i]+probabilities[i-1];
    }


}


int Protein::count_contacts(){
    //считается число топологических контактов НН
    int hh = 0;
    int position;
    for (int i =1; i<sequence.size()-1; i++){

        for ( std::pair <int, int> step : map_of_contacts[map_coordinate_to_int[conformation[i]]] ){


            if ( step!=conformation[i-1] && step!=conformation[i+1] && std::find(conformation.begin(), conformation.end(), step) !=conformation.end()  ){

                position=std::distance(conformation.begin(),find(conformation.begin(), conformation.end(), step));
                hh=hh+sequence[i]*sequence[position];
            }



        }


    }

    //концы конформации обрабатываются отдельно
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


int Protein::distance( std:: pair <int, int> point1, std:: pair <int, int>point2   ){
    //функция для подсчета расстояния между двумя точками на плоскости с учетом граничных условий
    int p1 = abs(point2.first-point1.first);
    int p2 = sequence.size()+3 -  abs(point2.first-point1.first);
    int part1 = std::min(p1, p2);
    p1 = abs(point2.second-point1.second);
    p2 = sequence.size()+3 -  abs(point2.second-point1.second);
    int part2 = std::min(p1, p2);

    return  part1 + part2;

}



void Protein::regrowth_end(int l ){
    //для выращивания участка длины l на конце цепочки

    std::vector<std::pair <int, int>>  C_t ;
    std::vector<int> seq_t;

    C_t.reserve(sequence.size());
    C_t.resize(sequence.size()-l);
    std::copy(conformation.begin(), conformation.begin()+sequence.size()-l, C_t.begin());


    seq_t.reserve(sequence.size());

    seq_t.resize(sequence.size()-l);
    std::copy(sequence.begin(), sequence.begin()+sequence.size()-l, seq_t.begin());


    int current_energy = count_contacts_breaked(seq_t, C_t, map_of_contacts, map_coordinate_to_int);
    int temp_e;
    static std::vector <std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0.0, 0));
    std::vector<int> energies;
    energies.resize(4, 0);
    float sum_probabilities=0.0;

    seq_t.emplace_back( sequence[sequence.size()-l]);


    for (int i =0; i<4; i++ ){


        //проверка на самопересечения + проверка первого шага, чтобы не начать строить в том же направлении, что и в старой конформации
        if (map_of_contacts[map_coordinate_to_int[C_t.back()]][i] == conformation[sequence.size()-l] ||  std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t.back()]][i]   )!=C_t.end()  ){
            //так как точка не подходит, то вероятность попасть в нее равна нулю
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        else{
            //точка подходит

            // Лучше потом переделать функцию для энергии
            C_t.emplace_back(map_of_contacts[map_coordinate_to_int[C_t.back()]][i]);
            temp_e = count_contacts_breaked(seq_t, C_t, map_of_contacts,map_coordinate_to_int);
            energies[i] = temp_e;
            probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

            // это все  energy-guided

            C_t.pop_back(); //????
        }
        sum_probabilities=sum_probabilities+probabilities_to_move[i].first;
    }


    if(sum_probabilities==0.0){
     // значит двигаться некуда, попытка перестройки оказалась неуспешной
        return;

    }
    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;
        // чтобы вероятности стали вероятностями

    }

    for (int i =1; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;
        // чтобы осуществить выбор с нужными вероятностями
    }

    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator); // случайно число от 0 до 1
    for (int i =0; i<4; i++ ){
        //забираем соответствующую координату
        if (q<probabilities_to_move[i].first){
            C_t.emplace_back( map_of_contacts[map_coordinate_to_int[ C_t.back()  ]][probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];

            break;

        }

    }

    // выстраиваем остальную часть цепочки
    for (int t = sequence.size()-l+1; t<sequence.size(); t++) {

        seq_t.emplace_back( sequence[t]);

        sum_probabilities = 0.0;


        for (int i =0; i<4; i++ ){
            // тут проверка только на самопересечения
            if(std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t.back()]][i])==C_t.end() ){

                C_t.emplace_back(map_of_contacts[map_coordinate_to_int[C_t.back()]][i]);
                temp_e = count_contacts_breaked(seq_t, C_t, map_of_contacts,map_coordinate_to_int);
                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

                C_t.pop_back();

            }
            else {
                //точка не подходит

                probabilities_to_move[i] = std::make_pair(0.0, i);
                energies[i] = 0;


            }

            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;

        }


        if(sum_probabilities==0.0){

            return;

        }

        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


        }
        for (int i =1; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;


        }

        q = distribution(generator);
        for (int i =0; i<4; i++ ){

            if (q<probabilities_to_move[i].first){
                C_t.emplace_back(map_of_contacts[map_coordinate_to_int[ C_t.back()  ]][probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];

                break;

            }

        }

    }


// дополнительная проверка, что раз мы дошли до этой части, участой нужной длины достроен
// по сути не нужно проверять, но мне так спокойнее
    if( C_t.size()==sequence.size()  ){


        q =distribution(generator);

        //это как раз вероятность, с которой будет приниматься новая конфромация
        // в оригинале это p=min( 1,   exp(..)  )
        // но так как q - число от 0 до 1, можно обойтись без функции min
        float probability_to_accept = exp(-(current_energy-E)/T);

        if(q<probability_to_accept){


            E = current_energy;
            conformation=C_t;


            if(E<min_E){
                results.clear();
                results.emplace_back(conformation);
                min_E=E;
                std:: cout << "New minimum fonded: " << min_E << std::endl ;

                iterator=iterator+1;
                if(change_T== false){
                    change_T=true;
                }
                std::ofstream out;
                out.open("coordinates_first_result.txt");

                for (auto c : conformation) {

                    out << c.first << " " << c.second << "   ";


                    out << std::endl;


                }
                out.close();




            }
            if (E==min_E){
                results.emplace_back(conformation);
            }


        }



    }
    else{
        //сюда не должны попадать никогда
        std::cout << " fail int the end " << std:: endl;
        std:: cout << "new size " << C_t.size() << std:: endl;
        return ;
    }


}



void Protein::regrowth_start(int l ) {

    std::vector<std::pair <int, int>>  C_t ;
    std::vector<int> seq_t;

    C_t.reserve(sequence.size());
    C_t.resize(sequence.size()-l);
    std::copy(conformation.begin()+l, conformation.end(), C_t.begin());



    seq_t.reserve(sequence.size());

    seq_t.resize(sequence.size()-l);
    std::copy(sequence.begin()+l, sequence.end(), seq_t.begin());


    int current_energy = count_contacts_breaked(seq_t, C_t, map_of_contacts, map_coordinate_to_int);
    int temp_e;
    static std::vector <std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0.0, 0));
    std::vector<int> energies;
    energies.resize(4, 0);
    std::pair<int, int> point;
    float sum_probabilities=0.0;

    seq_t.insert(seq_t.begin(), sequence[l-1]);

    for (int i =0; i<4; i++ ){
        //убираю вариант старой точки и самопересечения
        if (map_of_contacts[map_coordinate_to_int[C_t[0]]][i] == conformation[l-1] ||  std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t[0]]][i]   )!=C_t.end()  ){
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        else{
            //точка подходит

            // Лучше потом переделать функцию для энергии
            C_t.insert(C_t.begin(),map_of_contacts[map_coordinate_to_int[C_t[0]]][i]  );

            temp_e = count_contacts_breaked(seq_t, C_t, map_of_contacts,map_coordinate_to_int);
            energies[i] = temp_e;
            probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

            C_t.erase(C_t.begin());



        }



        sum_probabilities=sum_probabilities+probabilities_to_move[i].first;




    }


    if(sum_probabilities==0.0){
        //неудачная попытка

        return;

    }

    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


    }

    for (int i =1; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;

    }

    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator);
    for (int i =0; i<4; i++ ){

        if (q<probabilities_to_move[i].first){
            C_t.insert(C_t.begin(),map_of_contacts[map_coordinate_to_int[ C_t[0]  ]][probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];


            break;

        }

    }

    for (int t = l-2; t>-1; t--){
        seq_t.insert(seq_t.begin(), sequence[t]);

        sum_probabilities = 0.0;
        for (int i =0; i<4; i++ ){


            if(std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t[0]]][i])==C_t.end() ){


                C_t.insert(C_t.begin(),map_of_contacts[map_coordinate_to_int[C_t[0]]][i]  );

                temp_e = count_contacts_breaked(seq_t, C_t, map_of_contacts,map_coordinate_to_int);
                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

                C_t.erase(C_t.begin());



            }
            else {
                //точка не подходит


                probabilities_to_move[i] = std::make_pair(0.0, i);
                energies[i] = 0;

            }

            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;


        }


        if(sum_probabilities==0.0){

            return;

        }

        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


        }
        for (int i =1; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;


        }

        q = distribution(generator);

        for (int i =0; i<4; i++ ){

            if (q<probabilities_to_move[i].first){
                C_t.insert(C_t.begin(),map_of_contacts[map_coordinate_to_int[ C_t[0]  ]][probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];


                break;

            }

        }



    }


    if( C_t.size()==sequence.size()  ){

        q =distribution(generator);
        float probability_to_accept = exp(-(current_energy-E)/T);

        if(q<probability_to_accept){

            E = current_energy;
            conformation=C_t;

            iterator=iterator+1;

            if(change_T== false){
                change_T=true;
            }
            if(E<min_E){
                results.clear();
                results.emplace_back(conformation);
                min_E=E;
                std:: cout << "New minimum fonded: " << min_E << std::endl ;


                std::ofstream out;
                out.open("coordinates_first_result.txt");

                for (auto c : conformation) {

                    out << c.first << " " << c.second << "   ";


                    out << std::endl;


                }
                out.close();

            }
            if (E==min_E){
                results.emplace_back(conformation);
            }

        }

    }
    else{
        std::cout << " fail int the end " << std:: endl;
        std:: cout << "new size " << C_t.size() << std:: endl;
        return ;
    }


}


void Protein::regrowth_middle(int l, int start_position){


    int end_position = start_position+l-1;
    std::vector<std::pair <int, int>>  C_t, C_t_temp ;
    std::vector<int> seq_t, seq_t_temp;
    //std::copy(sequence.begin(), seq_t.begin()+ start_position, C_t_temp.begin());
    C_t.reserve(sequence.size());
    C_t_temp.resize(start_position);

    std::copy(conformation.begin(), conformation.begin()+ start_position, C_t_temp.begin());
    C_t = C_t_temp;


    C_t_temp.resize(sequence.size()-start_position-l);
    std::copy(conformation.begin()+end_position+1, conformation.end(), C_t_temp.begin());

    C_t.insert(C_t.end(), C_t_temp.begin(), C_t_temp.end());
    C_t_temp.clear();
    //это удаление куска конформации с start_position до end_position


    seq_t.reserve(sequence.size());
    seq_t_temp.resize(start_position);
    std::copy(sequence.begin(), sequence.begin()+ start_position, seq_t_temp.begin());

    seq_t = seq_t_temp;

    seq_t_temp.resize(sequence.size()-start_position-l);
    //seq_t_temp.clear();
    std::copy(sequence.begin()+end_position+1, sequence.end(), seq_t_temp.begin());
   // seq_t.reserve(seq_t.size()+seq_t_temp.size() );
    seq_t.insert(seq_t.end(), seq_t_temp.begin(), seq_t_temp.end());
    seq_t_temp.clear();
    //это удаление куска последовательности



    int current_energy = count_contacts_breaked(seq_t, C_t, map_of_contacts, map_coordinate_to_int);

    seq_t.insert(seq_t.begin()+start_position, sequence[start_position]);


    int temp_e;
    static std::vector <std::pair<float, int>> probabilities_to_move;
    probabilities_to_move.resize(4, std::make_pair(0.0, 0));


    std::vector<int> energies;
    energies.resize(4, 0);
    std::pair<int, int> point;
    float sum_probabilities=0.0;

    for (int i =0; i<4; i++ ){

         if (map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i] == conformation[start_position]){
            probabilities_to_move[i] = std::make_pair(0.0, i);
            energies[i] = 0;//strange, but it for time economy
            continue;
        }
        else  if( std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i]   )==C_t.end() &&distance(map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i], conformation[end_position+1]) <=abs(end_position+1-start_position)  ){


            // Лучше потом переделать функцию для энергии
            C_t.insert(C_t.begin()+start_position,map_of_contacts[map_coordinate_to_int[C_t[start_position-1]]][i]  );

            //temp_e = count_contacts_breaked(seq_t, C_t, map_of_contacts,map_coordinate_to_int);
            energies[i] =  temp_e = count_contacts_partied_t(seq_t, C_t, map_of_contacts,map_coordinate_to_int, start_position, current_energy);

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
        return;

    }

    sort(probabilities_to_move.begin(), probabilities_to_move.end());
    for (int i =0; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


    }

    for (int i =1; i<4; i++ ){
        probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;


    }


    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double q = distribution(generator);

    for (int i =0; i<4; i++ ){

        if (q<probabilities_to_move[i].first){
            C_t.insert(C_t.begin()+start_position,map_of_contacts[map_coordinate_to_int[ C_t[start_position-1]  ]][probabilities_to_move[i].second]  );
            current_energy = energies[probabilities_to_move[i].second];

            break;

        }

    }


    for (int t = start_position+1; t< end_position+1; t++){
       // std::cout << " t = " << t << " size = " << C_t.size() << std::endl;
        seq_t.insert(seq_t.begin()+t, sequence[t]);

        sum_probabilities = 0.0;


        for (int i =0; i<4; i++ ){

            if(std::find(C_t.begin(), C_t.end(), map_of_contacts[map_coordinate_to_int[C_t[t-1]]][i])==C_t.end() &&  distance( map_of_contacts[map_coordinate_to_int[C_t[t-1]]][i] ,conformation[end_position+1] )<=abs(end_position+1-t)){


                C_t.insert(C_t.begin()+t,map_of_contacts[map_coordinate_to_int[C_t[t-1]]][i]  );

                temp_e = count_contacts_partied_t(seq_t, C_t, map_of_contacts,map_coordinate_to_int, t, current_energy);

                energies[i] = temp_e;
                probabilities_to_move[i] = std::make_pair( exp(-(temp_e-current_energy)/T), i);

                C_t.erase(C_t.begin()+t);

            }
            else {
                //точка подходит

                probabilities_to_move[i] = std::make_pair(0.0, i);
                energies[i] = 0;


            }

            sum_probabilities=sum_probabilities+probabilities_to_move[i].first;

        }



        if(sum_probabilities==0.0){
            return;

        }

        sort(probabilities_to_move.begin(), probabilities_to_move.end());
        for (int i =0; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first/sum_probabilities;


        }
        for (int i =1; i<4; i++ ){
            probabilities_to_move[i].first= probabilities_to_move[i].first + probabilities_to_move[i-1].first;


        }




        q = distribution(generator);

        for (int i =0; i<4; i++ ){

            if (q<probabilities_to_move[i].first){
                C_t.insert(C_t.begin()+t,map_of_contacts[map_coordinate_to_int[ C_t[t-1]  ]][probabilities_to_move[i].second]  );
                current_energy = energies[probabilities_to_move[i].second];


                break;

            }

        }



    }

// проверка по сути не нужна, но мне так спокойнее

    if( C_t.size()==sequence.size()  ){

        q =distribution(generator);
        float probability_to_accept = exp(-(current_energy-E)/T);

        if(q<probability_to_accept){


            E = current_energy;
            conformation=C_t;

            if(E<min_E){
                results.clear();
                results.emplace_back(conformation);
                min_E=E;
                iterator=iterator+1;
                if(change_T== false){
                    change_T=true;
                }
                std:: cout << "New minimum fonded: " << min_E << std::endl ;

                std::ofstream out;
                out.open("coordinates_first_result.txt");

                for (auto c : conformation) {

                    out << c.first << " " << c.second << "   ";


                    out << std::endl;


                }
                out.close();




            }
            if (E==min_E){
                results.emplace_back(conformation);
            }





        }

    }
    else{
    std::cout << " fail int the end " << std:: endl;
    std:: cout << "new size " << C_t.size() << std:: endl;
    return ;
    }

}







void Protein::find_minimum() {

    double q;
    int l, start_position;
    std::default_random_engine generator(std::random_device{}() );
    std::uniform_real_distribution<double> distribution(0.0,1.0);


    std::default_random_engine generators1(std::random_device{}() );


    for (int iteration =1; iteration<1000000000000; iteration++) {


        q = distribution(generator);
        //std::cout<< q << " " ;

        for (int i = 0; i < probabilities.size(); i++) {
            if (q < probabilities[i]) {

                l =i+2;
                break;
            }


        }

        std::uniform_int_distribution<int> distribution1(0,sequence.size()-l);



        start_position = distribution1(generators1);

        if(start_position==0){
            regrowth_start(l);
        }
        else if (start_position==sequence.size()-l){
            regrowth_end(l);
        }
        else{
            regrowth_middle(l, start_position);
        }


        //std:: cout << l << " " ;

        if( iterator%50000==0 && change_T && T>0.1){
            T=T*0.98;
            std:: cout << "50000 moves are made; new T : " << T << std:: endl;
            change_T = false;
        }
        if(iteration%10000==0){
            std::cout << "Number of attempts :  " << iteration<< std::endl;
        }
        if(min_E==-9){
            std::cout << "-36 achieved! " << std:: endl;
            break;
        }



    }
    std::cout << std::endl;

};


int Protein::count_contacts_partied_t( std::vector <int> &sequence, std::vector <std::pair <int, int>> &conformation, std:: map <int, std::vector < std::pair <int, int> >> &map_of_contacts, std:: map <std::pair<int, int>, int> &map_coordinate_to_int, int t, int current_energy    ){


    int new_energy = current_energy;
    int position;
    for (std::pair<int, int> step : map_of_contacts[map_coordinate_to_int[conformation[t]]]){

        if (step != conformation[t - 1] && step != conformation[t + 1] &&
            std::find(conformation.begin(), conformation.end(), step) != conformation.end()) {

            position = std::distance(conformation.begin(), find(conformation.begin(), conformation.end(), step));
            new_energy = new_energy + sequence[t] * sequence[position];
        }




    }





    return new_energy ;

}




int count_contacts_breaked(std::vector <int> &sequence, std::vector <std::pair <int, int>> &conformation, std:: map <int, std::vector < std::pair <int, int> >> &map_of_contacts, std:: map <std::pair<int, int>, int> &map_coordinate_to_int    ) {
    int hh = 0;
    int position;
//static std::valarray<std::pair <int, int>>  steps = { std::make_pair(1, 0), std::make_pair(-1, 0), std::make_pair(0, 1),  std::make_pair(0, -1) };
//std:: vector <std::pair <int, int>> not_topological = {};
    std::pair<int, int> new_point, new_point_begin, new_point_end;
    for (int i = 1; i < sequence.size()-1; i++) {
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