/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Jiahao Miao
_/    _/  _/_/_/  _/_/_/_/ Matricola: 986136
*****************************************************************
*****************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>

#include "RndGen/random.h"

double error(double av, double av2, int n) {
  return n == 0 ? 0 : sqrt((av2 - av * av) / n);
}
double norm_squared(double x, double y, double z) {
  return x * x + y * y + z * z;
}

int main() {
  Random rnd;
  rnd.initialize();

  int M{100000};  // number of throws
  int N{100};     // number of blocks
  int L{M / N};   // throws in each bl ock
  int n{100};     // number of max steps
  int a{1};       // lattice constant

  double sum{};        // sum of the throws in each block
  double ave{};        // average of the throws in each block
  double av2{};        // average of the squares of the throws in each block
  double sum_prog{};   // sum of the averages
  double sum_prog2{};  // sum of the squares of the averages

  // double x{},y{},z{};
  std::ofstream fout{"data/3D_RW_discrete.dat"};
  double x[6] = {0., 0., 0.,
                 0., 0., 0.};  // position of the walker x,y,z,-x,-y,-z

  for (int i{1}; i <= n; i++) {  // max i steps
    sum_prog = 0;
    sum_prog2 = 0;
    for (int j{0}; j < N; j++) {  // data blocking
      sum = 0;
      for (int k{0}; k < L; k++) {  // L throws in each block
        x[0] = x[1] = x[2] = x[3] = x[4] = x[5] =
            0.;                       // initialize the position of the walker
        for (int l{0}; l < i; l++) {  // i steps to walk
          x[int(rnd.Rannyu(0, 6))] += a;  // randomly choose the
                                          // direction to walk
        }                                 // end of one walk
        /// position of the walker x,y,z,-x,-y,-z at the end of the walk
        /// from the origin
        sum += norm_squared(x[0] - x[3], x[1] - x[4], x[2] - x[5]);
      }
      ave = sqrt(sum / L);  // root mean square of the position
      av2 = ave * ave;      // squares of the averages
      sum_prog += ave;      // sum of the averages
      sum_prog2 += av2;     // sum of the squares of the averages
    }
    // print the results
    sum_prog /= double(N);   // average of the averages
    sum_prog2 /= double(N);  // average of the squares of the averages
    fout << i << " " << sum_prog << " " << error(sum_prog, sum_prog2, N - 1)
         << '\n';
  }

  fout.close();
  return 0;
}