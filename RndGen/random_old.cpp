/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random_old.h"

using namespace std;


void Random_old :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random_old :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random_old :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max)
   return min+(max-min)*Rannyu();
}

double Random_old :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random_old :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double Random_old :: Exponential(double lambda){
  // This function generates a random number from an exponential distribution with given lambda
  double s=Rannyu();
  return -log(1-s)/lambda;
}

double Random_old :: Lorentz(double mean, double gamma){
  // This function generates a random number from a Lorentzian distribution with given mean and gamma
  double s=Rannyu();
  return mean + gamma*tan(M_PI*(s-0.5));
}

void Random_old :: initialize(){
  // This function initializes the random number generator
  int seed[4]; // Array per memorizzare il seme
  int p1, p2; // Numeri primi per l'inizializzazione del generatore di numeri casuali

  // Lettura dei numeri primi da file
  std::ifstream Primes("Primes");
  if(Primes.is_open()){
    Primes >> p1 >> p2;
  }else{
    std::cerr << "PROBLEM: Impossibile aprire il file Primes" << std::endl;
    exit(1); // Uscita con errore
  }
  Primes.close();

  // Lettura del seme da file
  std::ifstream input("seed.in");
  std::string property;
  if (input.is_open()){
    while (!input.eof()){
      input >> property;
      if(property == "RANDOMSEED"){
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        SetRandom(seed, p1, p2);
      }
    }
    input.close();
  }else {
    std::cerr << "PROBLEM: Impossibile aprire il file seed.in" << std::endl;
    exit(1); // Uscita con errore
  }
}

double Random_old :: error(double av, double av2, int n){
  // This function calculates the standard deviation of the mean
  return n==0 ? 0 : sqrt((av2 - av*av)/n);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
