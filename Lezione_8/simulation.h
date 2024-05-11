#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <array>
#include <fstream>

#include "RndGen/random.h"

//Gaussian 
inline auto gauss = [](double x, double mu, double sigma){
    double xmu{(x - mu)/sigma};
    return exp(-xmu*xmu/2);
};

//Trial Wave Function
inline auto psi_T = [](double x, double mu, double sigma){
    return  gauss(x,mu,sigma) + gauss(x,-mu,sigma);
};

// Modulus squared of the trial wave function
inline auto psi_T2 = [](double x, double mu, double sigma){
    double psi{psi_T(x,mu,sigma)};
    return psi*psi; 
};

// Kinetic divided by the trial wave function 
inline auto Kin =  [](double x, double mu, double sigma){
    double xmu{x*mu};
    double sigma2{sigma*sigma};
    return -0.5/(sigma2*sigma2)* (x*x + mu*mu -sigma2 - 2*xmu*tanh( xmu /sigma2 ) );
};

// Potential Energy
inline auto Pot = [](double x){
    double x2{x*x};
    return (x2*x2 - 2.5*x2);
};

// Hamiltonian divided by the trial wave function
inline auto Hamiltonian = [](double x, double mu, double sigma){
    return Kin(x,mu,sigma) + Pot(x);
};

class Simulation{
    private:
        Random rnd; //random number generator
        int _naccept{},_nattempts{}; //for acceptance rate
        double _x, _delta,_delta_mu,_delta_sigma,_mu,_sigma ; //starting point,step size and parameters of the trial wave function 
        double _mu_error{},_sigma_error{}; //parameters of the trial wave function
        double _nblk{},_nsteps{},_steps{},_block_ave{},_block_err{};  //blocks, total steps, steps in one block, block average and block error
        double _Lnew{},_Lold{},_Lnew_err{},_Lold_err{}; //cost function for simulated annealing
        double _T{}; //temperature for simulated annealing
    public:
    
    //Constructor: <starting point, step size, mu, sigma> and initialize the random number generator
    Simulation();
    ~Simulation()=default;

    //input
    void input();

    //data blocking algorithm
    void data_blocking(int blocks, int steps);
    
    //Metropolis algorithm: returns the new value of x 
    double metropolis(std::function<double (double,double,double)> pdf, double x0, double delta);

    //Simulated Annealing algorithm: return the value of the cost function given a pdf and a cost function 
    double simulated_annealing(std::function<double (double,double,double)> pdf,std::function<double (double,double,double)> L, double T);

    //-integrate a function using metropolis algorithm to sample the pdf associated with the integral and return integral with nsteps montecarlo steps
    double integrate(std::function<double(double,double,double)> pdf,std::function<double(double,double,double)> f,int nsteps);
    
    //error function 
    inline double error(double ave, double ave2, int n) const{
        return n == 1 ? 0 : sqrt((ave2 - ave*ave)/(n-1));
    }
    //getters
    //return the value of mu
    inline double getmu() const {return _mu;}
    //return the value of sigma
    inline double getsigma() const {return _sigma;}
    //return the error of mu
    inline double getmu_error() const {return _mu_error;}
    //return the error of sigma
    inline double getsigma_error() const {return _sigma_error;}
    //return the block average
    inline double getblock_ave() const {return _block_ave;}
    //return the block error
    inline double getblock_err() const {return _block_err;}
    //return the acceptance rate
    inline double getacceptance() const {return static_cast<double>(_naccept)/_nattempts;}
    //return the value of x
    inline double getx() const {return _x;}
    //return the new value of the cost function
    inline double getLnew() const {return _Lnew;}
    //return the old value of the cost function
    inline double getLold() const {return _Lold;}
    //return the new value of the cost function error
    inline double getLnew_err() const {return _Lnew_err;}
    //return the old value of the cost function error
    inline double getLold_err() const {return _Lold_err;}
    //return the value of T
    inline double getT() const {return _T;}
    //return the value of delta_mu
    inline double getdelta_mu() const {return _delta_mu;}
    //return the value of delta_sigma
    inline double getdelta_sigma() const {return _delta_sigma;}
    //return the value of desired error for mu
    inline double getmu_error_desired() const {return _mu_error;}
    //return the value of desired error for sigma
    inline double getsigma_error_desired() const {return _sigma_error;}

};