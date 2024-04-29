/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Jiahao Miao
_/    _/  _/_/_/  _/_/_/_/ Matricola: 986136
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "RndGen/random.h"

int main(){
    Random rnd; // Oggetto per la generazione di numeri casuali
    rnd.initialize(); // Inizializzazione del generatore di numeri casuali (seed e numeri primi)

    // Apertura del file di output per scrivere i risultati
    std::ofstream fout;

    int N[4] = {1, 2, 10, 100}; // S_N 
    int M{10000}; // Numero di estrazioni casuali
    double sum{};
    for(int i = 0; i < 4; i++){
        fout.open("data/Lorentz_dice_" + std::to_string(N[i]) + ".dat");
        for(int j = 0; j < M; j++){
            sum = 0;
            for(int k = 0; k < N[i]; k++){
                sum += rnd.Lorentz(0, 1);
            }
            fout << sum / N[i] << std::endl;
        }
        fout.close();
    }
    
    return 0;
}

/**********************************************
*                 Parametri                   *
***********************************************
*     mu = 0                                  *
*     Gamma =1                                *
**********************************************/