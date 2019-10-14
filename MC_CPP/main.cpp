#include <iostream>
#include <fstream>
#include "mcmc.h"




int main() {
    //std::cout << "Hello, World!" << std::endl;
    std::ofstream out1;          // поток для записи
    out1.open("something.txt"); // окрываем файл для записи

    std::vector<float> result = {};
    std::vector<int> s_template = {1,1,1};


    float h =0.05;

    float start_t = 0.05;
    float current_t = start_t;
    float finish_t = 2.26;

    std::ofstream out2;          // поток для записи
    out2.open("conformations.txt");


    while(current_t<finish_t){
        Protein p(s_template);
        float f1 = p.MC_for_E(current_t, 0.55, 400000);

        result.push_back(f1);
        //current_t = current_t+h;
        out1<< f1<< " ";
        std::cout << current_t<< std::endl;

        current_t = current_t+h;

        for (int i = 0; i<p.conformation.size(); i++){
            out2 << p.conformation[i] << " ";


        }
        out2 << std::endl;
    }

    out1.close();
    //std::cout << f1<< std:: endl;

    return 0;


}