#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <array>

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
        double _x, _delta,_mu,_sigma ; //starting point,step size and parameters of the trial wave function 
        double _mu_error{},_sigma_error{}; //parameters of the trial wave function
        double _block_ave{},_block_err{};  //for data blocking
        double _Lnew{},_Lold{}; //cost function for simulated annealing
        double _montecarlo_error{}; //error of the montecarlo integration
    public:
    
    //Constructor: <starting point, step size, mu, sigma> and initialize the random number generator
    Simulation(double x0,double delta,double mu, double sigma);
    ~Simulation()=default;

    //data blocking algorithm
    void data_blocking(int blocks, int steps);
    //error function 
    inline double error(double ave, double ave2, int n)const;
    
    //Metropolis algorithm: returns the new value of x 
    double metropolis(std::function<double (double,double,double)> pdf, double x0, double delta);

    //Simulated Annealing algorithm: return the value of the cost function given a pdf and a cost function 
    double simulated_annealing(std::function<double (double,double,double)> pdf,std::function<double (double,double,double)> L, double T, int steps,std::array<double,2> delta);

    //-integrate a function using metropolis algorithm to sample the pdf associated with the integral and return integral with nsteps montecarlo steps
    //-it also saves the error of the montecarlo integration
    double integrate(std::function<double(double,double,double)> pdf,std::function<double(double,double,double)> f,int nsteps);
    
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
    //return the error of the montecarlo integration
    //the error is calculated by using the central limit theorem
    inline double getmontecarlo_error() const {return _montecarlo_error;}


};