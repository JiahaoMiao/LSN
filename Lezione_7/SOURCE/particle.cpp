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
#include <math.h>
#include "particle.h"

using namespace std;

void Particle :: initialize(){
   _spin = 1;
   _x.resize(_ndim);
   _xold.resize(_ndim);
   _v.resize(_ndim);
   return;
}

void Particle :: translate(vec delta, vec side){
   for(unsigned int i=0; i<_ndim; i++){
     _x(i) = pbc(_x(i) + delta(i), side(i));
   }
}

void Particle :: flip(){
   _spin = -1*this->getspin();
}

void Particle :: moveback(){
   _x = _xold;
}

void Particle :: acceptmove(){
   _xold = _x;
}

int Particle :: getspin(){
   return _spin;
}

void Particle :: setspin(int spin){
   _spin = spin;
   return;
}

double Particle :: getposition(int dim, bool xnew){
   if(xnew) return _x(dim);
   else return _xold(dim);
}

void Particle :: setposition(int dim, double position){
   _x(dim) = position;
   return;
}

void Particle :: setpositold(int dim, double position){
   _xold(dim) = position;
   return;
}

double Particle :: getvelocity(int dim){
   return _v(dim);
}

vec Particle :: getvelocity(){
   return _v;
}

void Particle :: setvelocity(int dim, double velocity){
   _v(dim) = velocity;
   return;
}

double Particle :: pbc(double position, double side){
  return position - side * rint(position / side);
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
