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

// Funzione per calcolare l'errore statistico
double error(double av, double av2, int n){
    return n == 0 ? 0 : sqrt((av2 - av*av)/n);
}

int main(){
    Random rnd; // Oggetto per la generazione di numeri casuali
    rnd.initialize(); // Inizializzazione del generatore di numeri casuali (seed e numeri primi)

    // Apertura del file di output per scrivere i risultati
    std::ofstream fout{"data/variance.dat"};

    // Parametri per il calcolo della varianza
    int M{100000}; // Numero totale di estrazioni casuali
    int N{100}; // Numero di blocchi
    int L{int(M / N)}; // Numero di estrazioni per blocco
    double ave = 0, av2 = 0, sum = 0, sum_prog = 0, sum2_prog = 0;

    // Calcolo della varianza
    for(int i = 0; i < N; i++){
        sum = 0;
        for(int j = 0; j < L; j++){
            sum += pow(rnd.Rannyu() - 0.5 , 2); // Aggiungo scarto quadratico
        }
        ave = sum / double(L); // Calcolo della media
        av2 = ave * ave; // Calcolo della media al quadrato
        sum_prog += ave; // Aggiunta della media alla somma progressiva
        sum2_prog += av2; // Aggiunta della media al quadrato alla somma progressiva dei quadrati

        // Scrittura dei risultati su file
        fout << (i + 1) * L << " " << sum_prog / double(i + 1) << " " << error(sum_prog / double(i + 1), sum2_prog / double(i + 1), i) << std::endl;
    }
    fout.close(); // Chiusura del file di output

    return 0;
}