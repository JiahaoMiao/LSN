#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "RndGen/random.h"
#include "TSP.h"


using namespace std;

int main(){

    vector<city> cities;
    ifstream Input("INPUT/cap_prov_ita.dat");
    if(!Input){
        cerr << "Error opening the input file" << "\n";
        return -1;
    }
    int N{};
    double x{},y{};
    while(Input >> x >>y){
        cities.emplace_back(city(x,y));
        N++;
    }
    cities.reserve(N);
    cities.shrink_to_fit();
    Input.close();

    // vector<int> way{19, 31, 77, 90, 23, 79, 55, 107, 24, 29, 27, 54, 92, 46, 17, 8,
    //    9, 75, 85, 7, 35, 20, 11, 59, 22, 42, 25, 70, 44, 37, 45, 83, 109,
    //    68, 94, 81, 93, 5, 32, 51, 2, 69, 82, 36, 78, 33, 84, 103, 100, 39,
    //    101, 74, 10, 99, 64, 108, 98, 15, 91, 16, 106, 52, 14, 57, 80, 66,
    //    28, 12, 47, 26, 102, 104, 60, 105, 13, 3, 95, 30, 41, 88, 6, 1, 38,
    //    67, 56, 58, 49, 71, 43, 53, 50, 72, 73, 76, 34, 4, 89, 40, 48, 62,
    //    87, 63, 86, 21, 18, 61, 96, 97, 65};

    vector<int> way{19, 31, 77, 90, 23, 79, 55, 107, 24, 29, 27, 54, 92, 46, 17, 8, 9, 75, 85, 7, 35, 20, 11, 59, 22, 42, 25, 70, 44, 37, 45, 83, 109, 68, 94, 81, 93, 5, 32, 51, 2, 69, 82, 36, 78, 33, 14, 57, 80, 66, 52, 106, 108, 64, 84, 103, 100, 39, 101, 74, 99, 10, 15, 98, 16, 91, 12, 47, 26, 102, 104, 60, 105, 13, 3, 95, 30, 41, 88, 6, 1, 38, 67, 56, 58, 49, 28, 71, 43, 53, 50, 72, 73, 76, 34, 4, 89, 40, 48, 62, 87, 63, 86, 21, 18, 61, 96, 97, 65};
    
    path route{way};
    std::cout << "Path L1 Length: " << route.LengthL1(cities) << "\n";
    std::cout << "Path L2 Length: " << route.LengthL2(cities) << "\n";

    return 0;
}