#include <iostream>
#include <string>
#include <cmath>
#include <array>
#include <functional>
#include <fstream>

#include "simulation.h"

using namespace std;

int main(int argc, char* argv[]){

    if(argc != 2){
        cerr << "Usage: " << argv[0] << "<cooling schedule>\n" ;
        cerr << "Available cooling schedules:\n";
        cerr << "0 = exponential, 1 = geometric, 2 = logarithmic, 3 = linear\n";
        return 1;
    }

    if(atoi(argv[1]) < 0 || atoi(argv[1]) > 3){
        cerr << "Invalid cooling schedule, digit an integer between 0 and 3\n";
        return 1;
    }

    int cooling = atoi(argv[1]);
    Simulation sim{};
    
    // Exponential Cooling
    auto exponentialCooling = [](double initialTemperature, double alpha, int iteration) {
        return initialTemperature * std::pow(alpha, iteration);
    };

    // Linear Cooling
    auto linearCooling= [](double initialTemperature, double beta, int iteration) {
        return initialTemperature - beta * iteration;
    };

    // Logarithmic Cooling
    auto logarithmicCooling = [](double initialTemperature,double alpha, int iteration) {
        return initialTemperature / std::log(alpha*iteration + M_E);
    };

    // Geometric Cooling
    auto geometricCooling =[](double initialTemperature, double alpha, int iteration) {
        return initialTemperature / (1 + alpha * iteration);
    };

    // Simulated Annealing
    array<string,4> cooling_schedules = {"exponential","geometric","logarithmic","linear",};

    fstream Tout{"DATA/SA/"+cooling_schedules[cooling]+ "/temperature.dat",ios::out};
    fstream Lout{"DATA/SA/"+cooling_schedules[cooling] + "/energy.dat",ios::out};
    fstream paramsout{"DATA/SA/"+cooling_schedules[cooling] + "/parameters.dat",ios::out};

    int k{};//SA iteration
    double T{sim.getT()};
    const double delta_mu{sim.getdelta_mu()};
    const double delta_sigma{sim.getdelta_sigma()};
    const double err_mu{sim.getmu_error()};
    const double err_sigma{sim.getsigma_error()};
    
    // //at the start Lold = 0 and Lold_err = 0
    // double diff{fabs(sim.getLnew()-sim.getLold())};
    // double err{sqrt(pow(sim.getLnew_err(),2) + pow(sim.getLold_err(),2))};

    //stopping criterion: when the step size in mu or sigma is smaller than the desired error then stop
    //This does not guarantee that the error is smaller than the desired error, it's the temperature that is too low
    while( 2*delta_mu *T > err_mu or 2*delta_sigma*T > err_sigma){ 
        switch(cooling){
            case 0:
                T = exponentialCooling(T,0.9999,k);
                break;
            case 1:
                T = geometricCooling(T,1e-4,k);
                break;
            case 2:
                T = logarithmicCooling(T,5e-4,k);
                break;
            case 3:
                T = linearCooling(T,1e-6,k);
                
                break;
            default:
                cerr << "Invalid cooling schedule\n";
                return 1;
        }

        Tout << k << " " << T << "\n";

        sim.simulated_annealing(psi_T2,Hamiltonian,T);
        
        Lout << k << " " << sim.getblock_ave() << " " << sim.getblock_err()<<  "\n";
        paramsout << sim.getmu() << " " << sim.getsigma() << "\n";
        
        if(k%10 == 0){
            cerr << "Completed iteration " << k << " with T = " << T << " and (mu,sigma) = (" << sim.getmu() << " , "<< sim.getsigma() << ")\n";
        }
        k++;
        // diff = fabs(sim.getLnew()-sim.getLold());
        // err = sqrt(pow(sim.getLnew_err(),2) + pow(sim.getLold_err(),2));
        // if(diff < err ){
        //     cerr << "Difference between the new and old value of the cost function = " << diff << " is lower than the error = " << err << "\n"; 
        // }
    }

    cerr << "Completed Simulated Annealing\n";
    cerr << "Final temperature = " << T << "\n";
    cerr << "Minimum found: H = " << sim.getL_min() << " with (mu,sigma) = (" << sim.getmu_min() << " , " << sim.getsigma_min() << ")\n"; 
    //after reaching the stopping criterion, data block to get values of mu and sigma with error
    //the temperature is the last value of T that makes the stopping criterion true
    double mu{},sigma{},H{};
    double sum_H{},sum2_H{};

    //before starting set the value of mu and sigma that make the cost function minimum (found during the simulated annealing) 
    sim.setmu(sim.getmu_min());
    sim.setsigma(sim.getsigma_min());

    ofstream wave{"DATA/Minimum/"+cooling_schedules[cooling]+"/wave.dat",ios::out};
    ofstream Hmin{"DATA/Minimum/"+cooling_schedules[cooling]+"/energy.dat",ios::out};
    ofstream paramsmin{"DATA/Minimum/"+cooling_schedules[cooling]+"/parameters.dat",ios::out};

    cerr << "starting to data block values of mu, sigma with error. Moving around the minimum using T = 0.001\n";
    T = 0.01; //so that delta_mu*T = 0.001, so this is the error 
    
    
    for(int j{}; j < 10000; j++){
        sim.simulated_annealing(psi_T2,Hamiltonian,T);
        mu = sim.getmu();
        sigma = sim.getsigma();
        //write the values of mu and sigma in a file
        paramsmin << j << " " << mu << " " << sigma << "\n";
    }

    cerr << "New values of mu and sigma found (mu,sigma) = (" << sim.getmu_min() << " , " << sim.getsigma_min() << ")\n"; 
    
    sim.setmu(sim.getmu_min());
    sim.setsigma(sim.getsigma_min());

    cerr << "Starting to sample the wave function with those new parameters\n";
    //Sampling of the modulus square of the wave function
    double x{};
    for (int i{1}; i <= 1000000; i++){
        x = sim.metropolis(psi_T2,x,3);
        wave  << x << "\n";
    }
    cerr << "Completed sampling of the wave function\n";
    cerr << "Starting to data block the value of the Energy\n";
    //Data blocking to get the value of the Energy
    for (int i{1}; i <= 100; i++){
        H = sim.integrate(psi_T2, Hamiltonian, 10000);
        sum_H += H;
        sum2_H += H*H;
         //write the values of the cost function in a file
        Hmin << i << " " << sum_H/i << " " << sim.error(sum_H/i,sum2_H/i,i) << "\n";
    }
    cerr << "Completed data blocking of the Energy\n";
   
    sum_H /= 100;
    sum2_H /= 100;

    cout << '\n'
            << ",=======================================," << '\n'
            << "|  H:          " << sum_H << " +/- " << sim.error(sum_H,sum2_H,100) << " |\n"
            << "'======================================='" << '\n' << '\n';
    
    return 0;
}