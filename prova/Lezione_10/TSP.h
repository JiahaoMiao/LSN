#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <unordered_set>
#include "RndGen/random.h"


class city{
    private:
        double _x,_y;
    public:
        city():_x(0),_y(0){}
        city(double x, double y):_x(x),_y(y){}
        inline double getX()const{return _x;}
        inline double getY()const{return _y;}
        inline void setX(double x){_x = x;}
        inline void setY(double y){_y = y;}
};

class path{
    private:
        std::vector<int> _path;
    public:
        path() = default;
        ~path() = default;
        explicit path(const path& p):_path(p._path){}
        path(path&& p) noexcept:_path(std::move(p._path)){}
        path(std::vector<int> path);

        path& operator=(path& p){_path = p._path; return *this;}
        path& operator=(std::vector<int> path){_path = path; return *this;}
        path& operator=(const path& p){if(this != &p){_path = p._path;} return *this;}
        path& operator=(path&& p)noexcept{
            if(this != &p){
                _path = std::move(p._path); 
            }
            return *this;
        }
        
        //check if all IDs are unique
        bool CheckPath(const std::vector<int>& vec);

        //Calculate the distance between two cities using the L1 norm
        inline double distanceL1(city c, city d);
        //Calculate the distance between two cities using the L2 norm
        inline double distanceL2(city c, city d);
        //Calculate the total distance of the path using the L1 norm
        double LengthL1(const std::vector<city>& cities); 
        //Calculate the total distance of the path using the L2 norm
        double LengthL2(const std::vector<city>& cities);
        //Print the path
        void PrintPath();
        //swap of two cities
        void swap(int i, int j){std::swap(_path[i],_path[j]);}
        //shift of m contiguous cities of n positions starting from the i-th city
        void shift(int i,int m, int n);
        //permutation among m contiguous cities starting from the i-th city
        void Permutation(int m,int i,int sep);
        //inversion of the order of m contiguous cities starting from the i-th city
        void inversion(int m,int i);

        //access to the path
        inline std::vector<int> getPath()const{return _path;}
        inline int operator[](int i)const{return _path[i];}
};

//Traveling Salesman Problem
class TSP{
    private:
        //vector of cities ordered from 0 to _cities.size()-1
        std::vector<city> _cities;
        std::vector<path> _paths;
        Random _rnd;
    public:
    TSP();
    TSP(std::vector<city> cities);
    TSP(std::vector<city> cities,std::vector<path> paths);
    TSP(std::vector<city> cities,std::vector<path> paths,Random rnd);

    //Mutations

    //pair permutation of cities
    path& PairPermutation(path& way);
    //shift of m contiguous cities of n positions
    path& Shift(path& way,int m, int n);
    //permutation among m contiguous cities 
    path& Permutation(path& way,int m);
    //inversion of the order of m contiguous cities
    path& Inversion(path& way,int m);

    //Selections
    int Selection();

    //Genetic Algorithm
    //New generation of paths, update the population, old population are deleted
    void NewGeneration();

    //Simulated Annealing
    void SimulatedAnnealing(const double& beta,path& way);


    //Mutations, p1= pair permutation, p2=shift, p3=permutation, p4=inversion
    path& Mutate(path& way,double p1,double p2,double p3,double p4);

    //Print the path with length
    void Print();
    void PrintBest();

    //Average length of best half of the population
    double AverageLength();
    //Best length of the current population
    double BestLength();
    //Print to file the Best path in Cartesian coordinates
    void BestPath(const std::string& filename);

    //Getters
    inline std::vector<city> getCities()const{return _cities;}
    inline std::vector<path> getPaths()const{return _paths;}
    //access to the path
    path& operator[](int i){
        // Check if index is within bounds
        if (i < 0 || (unsigned)i >= _paths.size()) {
            throw std::out_of_range("Index out of range");
        }
        return _paths[i];
    }

};

