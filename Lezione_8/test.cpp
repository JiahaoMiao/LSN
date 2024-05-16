#include <iostream>
#include "simulation.h"

using namespace std;


int main(){

    Simulation sim{};
    sim.setmu(0.86508);
    sim.setsigma(0.627365);
    //cout << Hamiltonian(1,0.86508,0.627365) << endl;
    cout << sim.integrate(psi_T2,Hamiltonian,100000) <<endl;
}