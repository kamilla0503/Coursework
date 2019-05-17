#include <iostream>
#include <fstream>
#include "fress.h"
using namespace std;

int main(int argc, char *argv[]) {
    //std::cout << "Hello, World!" << std::endl;
    cout.precision(3);
    //valarray <int> sequence = { 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0};
    //vector<int> sequence = {0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0};

   // vector<int> sequence {0, 0, 1};
  // vector <int> sequence = {0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0};

  //  vector <int> sequence = { 0,1, 0,1,1,1,1,1,1, 0};
 //vector<int> sequence = {1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1};

  string filename, hp_sequence;

    filename=argv[1];
    std::ifstream in(filename);
    vector <int> sequence;
    getline(in, hp_sequence);
    for (char c: hp_sequence){
        if(c=='H'){
            sequence.push_back(1);
        }
        if(c=='P'){
            sequence.push_back(0);
        }

    }
    cout << "length is " << sequence.size() << endl;

//vector <int> sequence = {0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0};
    Protein p(sequence);

    /**for (int i =0; i<p.probabilities.size(); i++){
        cout << p.probabilities[i] << " ";
    }**/
    int e = p.count_contacts();
    cout << "Start energy : " << e << endl;
    //p.find_minimum();

   // int d = p.distance(make_pair(0, 0), make_pair(0, 62));
  //  cout << "distance = " << d << endl;

    //cout << e << endl;


  // p.regrowth_middle(12, 20);
  // p.regrowth_end(10);

    time_t start, end;
    time (&start);


    p.find_minimum();
    time (&end);

    double dif = difftime(end, start);
    printf ("Time = %lf \n", dif);
    cout <<  " Number of moves : " <<p.iterator << std:: endl;
    cout << "Founden minimum : " << p.min_E << endl;

    cout << endl;
    return 0;
}