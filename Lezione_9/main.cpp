#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "RndGen/random.h"
#include "TSP.h"


using namespace std;

int main(int argc, char* argv[]){

    if(argc != 2){
        cerr << "Usage: " << argv[0] << "0/1\n";
        cerr << "0: Circumference\n1: Square\n";
        return -1;
    }
    int choice = atoi(argv[1]);
    std::vector<city> cities;
    Random rnd;
    rnd.initialize();
    //parameters
    int N = 34;//number of cities
    int L = 3000;//number of paths
    int M =300;//number of generations
    string filename{};
    switch (choice){
    case 0:{
        for(int i = 0; i < N; i++){
            double theta = rnd.Rannyu(0,2*M_PI);
            cities.emplace_back(city{cos(theta),sin(theta)});
        }
        filename = "DATA/Circumference/";
        break;
        }
    case 1:
        for(int i = 0; i < N; i++){
            cities.emplace_back(city{rnd.Rannyu(-1,1),rnd.Rannyu(-1,1)});
        }
        filename = "DATA/Square/";
        break;
    default:
        cerr << "Invalid choice\n";
        cerr << "Usage: " << argv[0] << "0/1\n";
        cerr << "0: Circumference\n1: Square\n";
        exit(-1);
    }

    //Print initial cities
    ofstream config(filename+"config.dat");
    for(auto& c: cities){
        config << c.getX() << " " << c.getY() << "\n";
    }
    //to close the path
    config << cities[0].getX() << " " << cities[0].getY() << "\n";
    config.close();
    //start with the ordered path
    //path always starts from 0 so i omit it and work with the rest of the cities
    vector<int> way(N-1);
    for(int i = 1; i < N; i++){
        way[i-1] = i;
    }
    path route{way};
    vector<path> paths(L);
    int i{},j{};
    for(int k{};k<L;k++){
        i = rnd.Rannyu(0,N-1);
        j = rnd.Rannyu(0,N-1);
        while(i==j){
            j=rnd.Rannyu(0,N-1);
        }
        // i % (N-1) and j % (N-1) to avoid the value N-1
        route.swap(i%(N-1),j%(N-1));//swap of two cities
        paths[k] = route;
    }
    TSP tsp{cities,paths};
     
    std::cout << "Initial best path" << "\n";
    tsp.PrintBest();

    ofstream ave{filename+"AverageLength.dat"};
    ofstream best{filename+"BestLength.dat"};

    ave << 0 << " " << tsp.AverageLength() << "\n";

    tsp.BestPath(filename+"/config/BestPath0.dat");

    for(int i{1};i<=M;i++){
        tsp.NewGeneration();
        std::cout << "Generation " << i << " ";
        tsp.PrintBest();
        ave << i << " " << tsp.AverageLength() << "\n";
        best << i << " " << tsp.BestLength() << "\n";
        if(i%2==0){
            tsp.BestPath(filename+"/config/BestPath"+to_string(i)+".dat");
        }
    }
    //print the best path to file
    tsp.BestPath(filename+"BestPath.dat");
    return 0;
}
