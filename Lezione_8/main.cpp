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
    int N{10000};
    int blk{100};

    sim.data_blocking(blk,N);

    cout << "Energy using data blocking: " << sim.getblock_ave()  << '\n';
    cout << "Error with data blocking: " << sim.getblock_err() << '\n';
    
    double energy = sim.integrate(psi_T2,Hamiltonian,10000);

    cout << "Energy using 1 calculation: " << energy  << '\n';
    cout << "Error : " << sim.getmontecarlo_error() << '\n';
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