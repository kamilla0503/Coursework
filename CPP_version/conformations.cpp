//
// Created by kamilla on 09.12.18.
//

#include<vector>
#include<map>
#include <tuple>
#include <iostream>
#include <algorithm>
using namespace std;

vector <int> vector_for_distance(vector<tuple<int, int>> saw ){
    vector <int> result;
    vector<vector<tuple<int, int>>> left;
    vector<tuple<int, int>> steps;
    int dx, dy;
    left.push_back({make_tuple(0,1), make_tuple(-1, 0)});
    left.push_back({make_tuple(1, 0), make_tuple(0, 1)});
    left.push_back({make_tuple(0, -1), make_tuple(1, 0)});
    left.push_back({make_tuple(-1, 0), make_tuple(0, -1)});
    for(int i=2; i<saw.size();i++){
        dx = get<0>(saw[i-1])- get<0>(saw[i-2]);
        dy =  get<1>(saw[i-1]) - get<1>(saw[i-2]);
        steps.push_back(make_tuple(dx, dy));
        dx = get<0>(saw[i])- get<0>(saw[i-1]);
        dy =  get<1>(saw[i]) - get<1>(saw[i-1]);
        steps.push_back(make_tuple(dx, dy));
        if(get<0>(saw[i])== get<0>(saw[i-1]) && get<0>(saw[i-1])==get<0>(saw[i-2]) ||get<1>(saw[i])== get<1>(saw[i-1]) && get<1>(saw[i-1])==get<1>(saw[i-2])){
            result.push_back(0);
        }
        else if(find(left.begin(), left.end(), steps)!=left.end() ){
            result.push_back(-1);
        }
        else{
            result.push_back(1);
        }
        steps.clear();
    }
        return result;
}

int distance_between_saws(vector<tuple<int, int>> saw1, vector<tuple<int, int>> saw2){
    int v1 = 0;
    int v2 = 0;
    vector <int>  s1 =  vector_for_distance(saw1);
    vector <int>  s2 = vector_for_distance(saw2);
    for (int i=0;i<s1.size(); i++){
        v1=v1+ abs( s1[i]-s2[i]);
        v2=v2+ abs( s1[i]+s2[i]);
    }

    if(v1<=v2){
        return v1;
    }
    else{
        return v2;
    }
}

vector<vector<tuple<int, int>>>  filter_conformations(vector<vector<tuple<int, int>>>  saws ){
    vector<vector<tuple<int, int>>> result;
    result.push_back(saws[0]);
    int k=0;
    for(int i=1; i<saws.size(); i++){
        k=0;
        if( i%100==0){

            cout << i << " "<< endl;
        }
        for ( vector<tuple<int, int>> conformation : result    ){
            if(  distance_between_saws(conformation, saws[i] )== 0){
                k=-1;
                break;
            }
        }
        if(k==-1){
            continue;
        }
        else{
            result.push_back(saws[i]);
        }
    }
    return result;
}

vector<vector<tuple<int, int>>> get_all_conformations(int length){
    static vector<tuple<int,int>> steps = {make_tuple(1, 0), make_tuple(-1, 0), make_tuple(0, 1),  make_tuple(0, -1)};
    //vector<tuple<int, int>> steps = {make_tuple(1, 0), make_tuple(-1, 0), make_tuple(0, 1),  make_tuple(0, -1)};
    vector<vector<tuple<int, int>>> result;
    vector<tuple<int, int>> temp;
    tuple<int, int> new_point;
    if(length==3){
        result.push_back({make_tuple(0, 0), make_tuple(1, 0), make_tuple(2, 0)});
        result.push_back({make_tuple(0, 0), make_tuple(1, 0), make_tuple(1, 1)});
        return result;
    }
    else{
        result=get_all_conformations(length-1);
        vector<vector<tuple<int, int>>>  new_conformations;
        for (int i=0; i<result.size(); i++){
            for (tuple<int, int> step : steps){
                //new_point = make_tuple(get<0>(result[i][-1])+get<0>(step), get<1>(result[i][-1])+get<1>(step) );
                new_point = make_tuple(get<0>(result[i].back())+get<0>(step), get<1>(result[i].back())+get<1>(step) );
                if(find(result[i].begin(), result[i].end(), new_point)!=result[i].end() ){
                    continue;
                }
                temp = result[i];
                temp.push_back(new_point);
                new_conformations.push_back(temp);
                temp.clear();

            }



        }

        return filter_conformations(new_conformations);

    }


}

