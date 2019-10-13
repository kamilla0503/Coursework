#include <iostream>
#include <fstream>
#include "mcmc.h"




int main() {
    //std::cout << "Hello, World!" << std::endl;
    std::ofstream out1;          // поток для записи
    out1.open("something.txt"); // окрываем файл для записи

    std::vector<float> result = {};
    std::vector<int> s_template = {1,1,1};


    float h =0.1;

    float start_t = 0.05;
    float current_t = start_t;
    float finish_t = 2.26;
    while(current_t<finish_t){
        Protein p(s_template);
        float f1 = p.MC_for_E(current_t, 0.55, 400000);

        result.push_back(f1);
        //current_t = current_t+h;
        out1<< f1<< " ";
        std::cout << current_t<< std::endl;

        current_t = current_t+h;
    }

    out1.close();
    //std::cout << f1<< std:: endl;

    return 0;
}