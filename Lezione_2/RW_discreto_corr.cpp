/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Jiahao Miao
_/    _/  _/_/_/  _/_/_/_/ Matricola: 986136
*****************************************************************
*****************************************************************/

#include <array>
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

  constexpr unsigned int M{100000};  // number of throws
  constexpr unsigned int N{100};     // number of blocks
  constexpr unsigned int L{M / N};   // throws in each bl ock
  constexpr unsigned int n{100};     // number of max steps
  double a{1.};                      // lattice constant

  double sum{};        // sum of the throws in each block
  double ave{};        // average of the throws in each block
  double av2{};        // average of the squares of the throws in each block
  double sum_prog{};   // sum of the averages
  double sum_prog2{};  // sum of the squares of the averages

  // double x{},y{},z{};
  std::ofstream fout{"data/3D_RW_discrete_corr.dat"};
  // L random walks' finale position after some steps
  // N blocks
  std::array<std::array<double, L>, N> x{};
  std::array<std::array<double, L>, N> y{};
  std::array<std::array<double, L>, N> z{};

  std::array<double, 6> dir{};

  int appo{};

  for (unsigned int i{1u}; i <= n; i++) {  // i-th steps
    sum_prog = 0;
    sum_prog2 = 0;
    // in each block j, restart to walk from the i-th step
    for (unsigned int j{0u}; j < N; j++) {  // data blocking
      sum = 0;
      // L random walks for each block
      for (unsigned int k{}; k < L; k++) {
        appo = int(rnd.Rannyu(0, 6));
        dir[appo] = a;  // randomly choose the direction
        x[j][k] += dir[0] - dir[3];
        y[j][k] += dir[1] - dir[4];
        z[j][k] += dir[2] - dir[5];
        sum += norm_squared(x[j][k], y[j][k], z[j][k]);
        dir[appo] = 0;  // reset
      }                 // end of the L RW with i steps
      // reset after each block
      ave = sqrt(sum / L);  // root mean square of the position
      av2 = ave * ave;      // squares of the averages
      sum_prog += ave;      // sum of the averages
      sum_prog2 += av2;     // sum of the squares of the averages
      // print the results
    }
    fout << i << " " << sum_prog / double(N) << " "
         << error(sum_prog / double(N), sum_prog2 / double(N), N - 1) << '\n';
  }
  fout.close();
  return 0;
}