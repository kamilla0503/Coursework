#include <iostream>
#include <time.h>
//#include "conformations.cpp"
#include "functions.h"

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


    vector<vector<tuple<int, int>>> test = get_all_conformations(10);

    time (&end);
    cout << endl;
    cout << "number of conformations " << test.size() << endl;

    double dif = difftime(end, start);
    printf ("Time = %lf \n", dif);

    std::cout << "Hello, World!" << std::endl;
    return 0;
}