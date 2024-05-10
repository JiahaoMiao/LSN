#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <functional>
#include <fstream>

#include "simulation.h"

using namespace std;

int main(){
    double mu{1},sigma{0.5};
    double x0{0};
    double delta{1};
    Simulation sim{x0,delta,mu,sigma};
    // int N{10000};
    int blk{100};
    
    ofstream out{"DATA/blk100steps100.dat"};

    double ave{},ave2{},appo{};

    for(int i{1}; i <= blk; i++){
        appo = sim.integrate(psi_T2,Hamiltonian,100);
        ave+= appo;
        ave2+= appo*appo;
        out << ave/i << " " << sim.error(ave/i,ave2/i,i) << '\n';
    }

    out.close();
    out.open("DATA/blk100steps1000.dat");
    ave = 0;
    ave2=0;
    for(int i{1}; i <= blk; i++){
        appo = sim.integrate(psi_T2,Hamiltonian,1000);
        ave+= appo;
        ave2+= appo*appo;
        out << ave/i << " " << sim.error(ave/i,ave2/i,i) << '\n';
    }
    ave = 0;
    ave2=0;
    out.close();
    out.open("DATA/blk100steps10000.dat");
    for(int i{1}; i <= blk; i++){
        appo = sim.integrate(psi_T2,Hamiltonian,10000);
        ave+= appo;
        ave2+= appo*appo;
        out << ave/i << " " << sim.error(ave/i,ave2/i,i) << '\n';
    }
    ave = 0;
    ave2=0;
    out.close();
    out.open("DATA/blk100steps100000.dat");
    for(int i{1}; i <= blk; i++){
        appo = sim.integrate(psi_T2,Hamiltonian,100000);
        ave+= appo;
        ave2+= appo*appo;
        out << ave/i << " " << sim.error(ave/i,ave2/i,i) << '\n';
    }

    out.close();

    // double T{10000};

    // int steps{1000};

    // // Exponential Cooling
    // auto exponentialCooling = [](double initialTemperature, double alpha, int iteration) {
    //     return initialTemperature * std::pow(alpha, iteration);
    // };

    // // Linear Cooling
    // auto linearCooling= [](double initialTemperature, double beta, int iteration) {
    //     return initialTemperature - beta * iteration;
    // };

    // // Logarithmic Cooling
    // auto logarithmicCooling = [](double initialTemperature, int iteration) {
    //     return initialTemperature / std::log(iteration + 1);
    // };

    // // Boltzmann Cooling
    // auto boltzmannCooling = [](double initialTemperature, int iteration) {
    //     return initialTemperature / std::log(iteration + std::exp(1));
    // };

    // // Geometric Cooling
    // auto geometricCooling =[](double initialTemperature, double alpha, int iteration) {
    //     return initialTemperature / (1 + alpha * iteration);
    // };

    // Simulated Annealing

    // Exponential Cooling




    return 0;
}