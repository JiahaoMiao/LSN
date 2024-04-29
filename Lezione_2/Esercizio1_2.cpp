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

double p_x(double x) {
  return 1 - sqrt(1 - x);  // inverse of the cumulative distribution function
                           // of 2(1-x) in [0,1]
}

int main() {
  Random rnd;
  rnd.initialize();

  // integrale di f(x) = pi/2 * cos(pi/2 * x) con x in [0,1]

  int M{100000};  // number of throws
  int N{100};     // number of blocks
  int L{M / N};   // throws in each block

  double sum{};        // sum of the throws in each block
  double ave{};        // average of the throws in each block
  double av2{};        // average of the squares of the throws in each block
  double sum_prog{};   // sum of the averages
  double sum_prog2{};  // sum of the squares of the averages

  std::ofstream fout{"data/importance_sample.dat"};

  for (int i = 0; i < N; i++) {
    sum = 0;
    for (int j{0}; j < L; j++) {
      double x = p_x(rnd.Rannyu());
      sum += M_PI / 2 * cos(M_PI / 2 * x) / (2 * (1 - x));
    }
    ave = sum / L;     // integral in each block
    av2 = ave * ave;   // squares of the averages
    sum_prog += ave;   // sum of the averages
    sum_prog2 += av2;  // sum of the squares of the averages

    fout << (i + 1) << " " << sum_prog / double(i + 1) << " "
         << error(sum_prog / double(i + 1), sum_prog2 / double(i + 1), i)
         << '\n';
  }

  return 0;
}