/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Jiahao Miao
_/    _/  _/_/_/  _/_/_/_/ Matricola: 986136
*****************************************************************
*****************************************************************/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include "RndGen/random.h"

int main(){
    Random rnd;
    rnd.initialize();

    int M {1000000};
    int N {100};
    int L{M/N};
    double ave {};
    double ave2 {};
   
    double C_prog {0};
    double C_prog2 {0};

    double P_prog {0};
    double P_prog2 {0};

    std::ofstream call{"data/Call_discrete.dat"};
    std::ofstream put{"data/Put_discrete.dat"};
    
    //Parameters for the call option
    double S0 {100};
    double T {1};
    double K {100};
    double r {0.1};
    double sigma {0.25};

    double Z{};
    double S{};
    double C {0};
    double P {0};
    
    for(int i{0}; i<N; i++){ //data blocking
        C = 0;
        P = 0;
        for(int j{0}; j<L; j++){
            S = S0;
            for(int k{0}; k<100; k++){ //GBM
                Z = rnd.Gauss(0,1);
                S *= exp((r-0.5*sigma*sigma)*T/100. + sigma*Z*sqrt(T/100.));
            } 
            C += exp(-r*T)*std::max(0.,S-K);
            P += exp(-r*T)*std::max(0.,K-S);
        }
        ave = C/L;
        ave2 = ave*ave;

        C_prog += ave;
        C_prog2 += ave2;

        ave = P/L;
        ave2 = ave*ave;

        P_prog += ave;
        P_prog2 += ave2;

        call << (i+1)*L << " " << C_prog/double(i+1) << " " << rnd.error(C_prog/double(i+1), C_prog2/double(i+1), i) << std::endl;
        put << (i+1)*L << " " << P_prog/double(i+1) << " " << rnd.error(P_prog/double(i+1), P_prog2/double(i+1), i) << std::endl;
    }

    call.close();
    put.close();
    return 0;
}