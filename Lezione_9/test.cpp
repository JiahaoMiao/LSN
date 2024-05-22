#include "TSP.h"

using namespace std;

int main(){

    std::vector<city> cities;
    Random rnd;
    rnd.initialize();
    double theta{2*M_PI/34};//angle between two cities
    for(int i = 0; i < 34; i++){
        cities.emplace_back(city{cos(theta),sin(theta)});
    }

    //start with the ordered path
    //path always starts from 0 so i omit it and work with the rest of the cities
    vector<int> way(33);
    for(int i = 1; i < 34; i++){
        way[i-1] = i;
    }
    path p{way};
    
    std::vector<path> paths;
    paths.emplace_back(p);

    TSP tsp{cities,paths};
    
    cout << "Swap of 2 random different cities \n";

    tsp.getPaths()[0].PrintPath();
    p = tsp.PairPermutation(0);
    p.PrintPath();

    cout << "Shift of 3 cities of 2 positions starting from the 4-th city\n";
    path p1{way};
    p1.PrintPath();
    p1.shift(4,3,2);//shift 3 cities of 2 positions starting from the 4-th city
    p1.PrintPath();

    cout << "Permutation of 4 cities with other 4 different cities starting from the 5-th city\n";

    path p2{way};
    p2.PrintPath();
    p2.Permutation(4,5,3);//Permutation of 4 cities with other 4 different cities starting from the 5-th city
    p2.PrintPath();

    cout << "Inversion of the order of 5 contiguous cities starting from the 6-th city\n";
    path p3{way};
    p3.PrintPath();
    p3.inversion(5,6);//Inversion of the order of 5 contiguous cities starting from the 6-th city
    p3.PrintPath();
    return 0;
}
