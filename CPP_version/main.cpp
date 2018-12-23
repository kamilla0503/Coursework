#include <iostream>
#include <time.h>
//#include "conformations.cpp"
#include "conformations.h"
//#include <fstream>
int main() {

    vector<tuple<int, int>> saw_axample, ssaw_axample ;
    saw_axample.push_back(make_tuple(0, 0));
    saw_axample.push_back(make_tuple(0, 1));
    saw_axample.push_back(make_tuple(0, 2));
    saw_axample.push_back(make_tuple(0, 3));
    saw_axample.push_back(make_tuple(1, 3));
    saw_axample.push_back(make_tuple(2, 3));
    saw_axample.push_back(make_tuple(2, 2));

    ssaw_axample.push_back(make_tuple(0, 0));
    ssaw_axample.push_back(make_tuple(0, -1));
    ssaw_axample.push_back(make_tuple(0, -2));
    ssaw_axample.push_back(make_tuple(0, -3));
    ssaw_axample.push_back(make_tuple(1, -3));
    ssaw_axample.push_back(make_tuple(2, -3));
    ssaw_axample.push_back(make_tuple(3, -3));


    vector <int> descritption = vector_for_distance(saw_axample);
    for (int i=0; i<descritption.size(); i++){

        cout << descritption[i]<< " ";

    }
    cout << endl;

    cout << "result of comparision " << distance_between_saws(saw_axample, ssaw_axample) << endl;

    time_t start, end;
    time (&start);

    int walk_len;
    std::cin >> walk_len;

    vector<vector<tuple<int, int>>> test = get_all_conformations(walk_len);

    time (&end);
    cout << endl;
    cout << "number of conformations " << test.size() << endl;


    //f=fopen("out_12.dat", "ab+");

    double dif = difftime(end, start);
    printf ("Time = %lf \n", dif);
    /**for (int i=0; i<test.size(); i++) {

        fwrite(test[i], sizeof(test[i]), 1, f);

    }**/


    /**
    ofstream fout("data.dat", ios::out | ios::binary);
    fout.write((char*)&student[0], student.size() * sizeof(Student));
    fout.close();
    **/

    std::cout << "Hello, World!" << std::endl;
    return 0;
}