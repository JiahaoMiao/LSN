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

    while( 2*delta_mu *T > err_mu or 2*delta_sigma*T > err_sigma){ // stopping criterion, when the difference between the new and old value of the cost function is less than the error
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
    //after reaching the stopping criterion, data block to get values of mu and sigma with error
    double mu{},sigma{},H{};
    double sum_mu{},sum2_mu{};
    double sum_sigma{},sum2_sigma{};
    double sum_H{},sum2_H{};
    for(int i{}; i < 100; i++){
        
        sim.simulated_annealing(psi_T2,Hamiltonian,T);
        mu = sim.getmu();
        sigma = sim.getsigma();
        H = sim.getblock_ave();
        sum_mu += mu;
        sum2_mu += mu*mu;
        sum_sigma += sigma;
        sum2_sigma += sigma*sigma;
        sum_H += H;
        sum2_H += H*H;
    }
    sum_mu /= 100;
    sum2_mu /= 100;
    sum_sigma /= 100;
    sum2_sigma /= 100;
    sum_H /= 100;
    sum2_H /= 100;

    cout << '\n'
            << ",=======================================," << '\n'
            << "|  H:          " << sum_H << " +/- " << sim.error(sum_H,sum2_H,100) << " |\n"
            << "| mu:          " << sum_mu <<   " +/- " << sim.error(sum_mu, sum2_mu, 100) << " |\n"
            << "| sigma:       " << sum_sigma  << " +/- " << sim.error(sum_sigma,sum2_sigma,100) << " |\n"
            << "'======================================='" << '\n' << '\n';
    
    return 0;
}