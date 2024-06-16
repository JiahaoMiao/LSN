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
#include "system.h"

using namespace std;

void displayProgressBar(int progress, int total, int width = 50) {
    float percentage = static_cast<float>(progress) / total;
    int barWidth = static_cast<int>(percentage * width);

    std::cout << "[";

    for (int i = 0; i < width; ++i) {
        if (i < barWidth)
            std::cout << "=";
        else
            std::cout << " ";
    }

    std::cout << "] " << std::fixed << std::setprecision(1) << (percentage * 100) << "%\r";
    std::cout.flush();
}

int main (int argc, char *argv[]){

  if(argc != 3){
    cerr << "Usage: " << argv[0] << " <phase> <Algorithm>\n";
    cerr << "<phase> :\n 0 = gas\n 1 = liquid\n 2 = solid\n 3 = Ising\n";
    cerr << "<Algorithm> :\n 0 = MD\n 1 = MC\n";
    return -1;
  }

  int phase{atoi(argv[1])};
  int algorithm{atoi(argv[2])};
  int steps_to_skip{};
  int nconf{1};
  System SYS;
  SYS.initialize(phase,algorithm);
  SYS.initialize_properties();
  SYS.block_reset(0);

  switch(phase){//to make sure the system is equilibrated
    case 0:
      steps_to_skip = 10000;
      break;
    case 1:
      steps_to_skip = 2000;
      break;
    case 2:
      steps_to_skip = 2000;
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
      displayProgressBar(i+1, SYS.get_nbl());
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
