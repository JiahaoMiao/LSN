/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Jiahao Miao
_/    _/  _/_/_/  _/_/_/_/ Matricola: 986136
*****************************************************************
*****************************************************************/
#include "wave_func.h"
 
int main() {
  WaveFunctionSimulation sim;
  sim.Input();

  // equilibrate
  for (int i = 0; i < 500; i++) sim.Move();

  for (int i = 0; i < sim.getnblk(); i++) { // data blocking
    sim.ResetAll();

    for (int j = 0; j < sim.getnstep(); j++) {// simulation
      sim.SavePos(); // current position
      sim.Move();
      sim.Accumulate();
    } 

    sim.PrintAccRate(i); // acceptance rate
    sim.BlockAverages();
    sim.SaveDist(i); // save out block-averaged distances
  }

  return 0;
}
