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

    return 0;
}