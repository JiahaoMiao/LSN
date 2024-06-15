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
#include <iomanip>
#include "system.h"

using namespace std;

void print_progress_bar(double progress) {
    int barWidth = 70;
    cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
}

int main (int argc, char *argv[]){

  if(argc != 2){
    cerr << "Usage: " << argv[0] << " <phase>\n";
    cerr << "<phase> :\n 0 = gas\n 1 = liquid\n 2 = solid\n 3 = Ising\n";
    return -1;
  }

  int phase{atoi(argv[1])};
  int steps_to_skip{};
  int nconf{1};
  System SYS;
  SYS.initialize(phase);
  SYS.initialize_properties();
  SYS.block_reset(0);

  switch(phase){//to make sure the system is equilibrated
    case 0:
      steps_to_skip = 15000;
      break;
    case 1:
      steps_to_skip = 400;
      break;
    case 2:
      steps_to_skip = 400;
      break;
    case 3:
      steps_to_skip = 4000;
      break;
    default:
      cerr << "Invalid phase\n";
      return -1;
  }

  for (int i{}; i < steps_to_skip; i++) SYS.step(); //Equilibration

  //data blocking
  for(int i{}; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j{}; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step(); //Move particle with verlet algorithm
      SYS.measure(); //Properties measurement
      // if(j%40 == 0){
      //   SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
      //   nconf++;
      // }
      }
      print_progress_bar((double)(i+1)/SYS.get_nbl());
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    SYS.finalize();

  return 0;
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
