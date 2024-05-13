#include "simulation.h"

Simulation::Simulation(){
    input();
    rnd.initialize();
    //initialize the cost function for simulated annealing
    data_blocking(_nblk,_steps);
    _Lnew = _block_ave;
    _Lnew_err = _block_err;
    //initialize the minimum value of the cost function
    _L_min = _Lnew;
    _mu_min = _mu;
    _sigma_min = _sigma;
}

void Simulation::input(){
    std::ifstream ReadInput;

    std::cerr << '\n'
        << ",=====================================," << '\n'
        << "| Simulated Annealing                 |" << '\n'
        << "| Monte Carlo simulation (Metropolis) |" << '\n'
        << "'====================================='" << "\n\n";
    
    //Read input informations
    ReadInput.open("input.dat");

    if(!ReadInput.is_open()){
        std::cerr << "Error: unable to open input.dat\n";
        std::exit(1);
    }

    ReadInput >> _nblk;
    ReadInput >> _nsteps;
    ReadInput >> _delta;
    ReadInput >> _x;

    ReadInput >> _mu;
    ReadInput >> _sigma;

    ReadInput >> _T; //initial temperature             

    ReadInput >> _delta_mu;     // step size for mu
    ReadInput >> _delta_sigma;    // step size for sigma
    ReadInput >> _mu_error; // desired error for mu
    ReadInput >> _sigma_error; // desired error for sigma

    _steps = _nsteps/_nblk;

    std::cerr << "Total number of steps = " << _nsteps << '\n';
    std::cerr << "Number of blocks = " << _nblk << '\n';
    std::cerr << "Number of steps in one block = " << _steps << '\n' << '\n';
    std::cerr << "Step lenght = " << _delta << '\n';
    std::cerr << "Initial position = " << _x << '\n';
    std::cerr << "Initial mu = " << _mu << '\n';
    std::cerr << "Initial sigma = " << _sigma << '\n';
    std::cerr << "Initial temperature = " << _T << '\n';
    std::cerr << "Step size for mu = " << _delta_mu << '\n';
    std::cerr << "Step size for sigma = " << _delta_sigma << '\n';
    std::cerr << "Desired error for mu = " << _mu_error << '\n';
    std::cerr << "Desired error for sigma = " << _sigma_error << '\n';

    ReadInput.close();
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
   
double Simulation::simulated_annealing(std::function<double (double,double,double)> pdf,std::function<double(double,double,double)> L,double T){
    
    _Lold = _Lnew;
    _Lold_err = _Lnew_err;
    
    _T = T;

    //save the old parameters in case the new ones are not accepted
    double mu_old{_mu};
    double sigma_old{_sigma};

    //change the parameters with uniform distribution
    //at low T the parameters change very little
    _mu += rnd.Rannyu(-_delta_mu,_delta_mu)*T;

    do{
    _sigma += rnd.Rannyu(-_delta_sigma,_delta_sigma)*T;
    if(fabs(_sigma) < 0.002) std::cerr << "Sigma "<< _sigma <<"\n";
    }while(fabs(_sigma) < 0.002); //sigma cannot be too small or the wave function diverges

    data_blocking(_nblk,_steps);
    _Lnew = _block_ave;
    _Lnew_err = _block_err;

    //save the minimum value of the cost function
    //if the new value is lower than the minimum value, then it becomes the new minimum so far
    if(_Lnew < _L_min){ 
        _L_min = _Lnew;
        _mu_min = _mu;
        _sigma_min = _sigma;
    }

    //Metropolis algorithm to accept or reject the new parameters
    double p = exp(-( _Lnew - _Lold )/T);

    if(p >= 1){
        return _Lnew;
    }else{
        if(rnd.Rannyu() < p){
            return _Lnew;
        }
    }
    //not accepted, return the old parameters
    _mu = mu_old;
    _sigma = sigma_old;
    return _Lold;
}