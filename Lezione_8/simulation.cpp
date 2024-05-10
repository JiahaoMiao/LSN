#include "simulation.h"

Simulation::Simulation(double x0,double delta,double mu, double sigma):_x{x0},_delta{delta},_mu{mu},_sigma{sigma}{
    rnd.initialize();
    _Lnew = integrate(psi_T2,Hamiltonian,10000); //initialize the cost function to prepare for simulated annealing
}

void Simulation::data_blocking(int nblocks,int nsteps){
    double ave{}, sum{}, sum2{};
    for(int i{}; i < nblocks; i++){
        ave = integrate(psi_T2, Hamiltonian, nsteps);
        sum += ave;
        sum2 += ave*ave;
    }
    _block_ave = sum/nblocks;
    _block_err = error(sum/nblocks,sum2/nblocks,nblocks);
}

double Simulation::metropolis(std::function<double (double,double,double)> pdf, double x0, double delta){
    //uniform random number generator
    _nattempts++;
    double x = x0 + rnd.Rannyu(-delta, delta);
    double p = pdf(x,_mu,_sigma) / pdf(x0,_mu,_sigma);
    if(p >= 1){
        _naccept++;
        return x;
    }else{
        if(rnd.Rannyu() < p){
            _naccept++;
            return x;
        }
        
    }
    return x0; //if the new value is not accepted, return the old value 
}

//can try inlining this function to see if it improves performance
double Simulation::integrate(std::function<double(double,double,double)> pdf,std::function<double(double,double,double)> f,int nsteps){
    double integral{};
    for(int i{}; i < nsteps; i++){
        //sample the modulus squared of the wave function
        _x = this->metropolis(pdf,_x,_delta); //does the same as calling directly metropolis(pdf,_x,_delta)
        integral += f(_x,_mu,_sigma); //integrate the function f    
    }
    
    return integral/nsteps;
}
   
double Simulation::simulated_annealing(std::function<double (double,double,double)> pdf,std::function<double(double,double,double)> L, double T, int steps,std::array<double,2> delta){
    
    _Lold = _Lnew;
    
    //randomly change the parameters of the trial wave function
    _mu += rnd.Rannyu(-delta[0],delta[0]);
    _sigma += rnd.Rannyu(-delta[1],delta[1]);

    _Lnew = integrate(pdf, L, steps);

    //Metropolis algorithm to accept or reject the new parameters
    double p = exp(-(_Lnew - _Lold)/T);

    if(p >= 1){
        return _Lnew;
    }else{
        if(rnd.Rannyu() < p){
            return _Lnew;
        }
    }
    return _Lold;
}