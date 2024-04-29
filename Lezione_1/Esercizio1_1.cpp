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
// Parametri:
// - av: media
// - av2: media al quadrato
// - n: numero di misure
// Restituisce l'errore statistico
double error(double av, double av2, int n){
    return n == 0 ? 0 : sqrt((av2 - av*av)/n);
}
 
int main(){
    Random rnd{}; // Inizializzazione dell'oggetto Random per generare numeri casuali
    rnd.initialize(); // Inizializzazione del generatore di numeri casuali (seed e numeri primi)
   
    // Parametri per il calcolo della media
    int M{100000}; // Numero totale di estrazioni casuali
    int N{100}; // Numero di blocchi
    int L{int( M / N)}; // Numero di estrazioni per blocco
    double ave{}, av2{}, sum{}, sum_prog{}, sum2_prog{};

    // Apertura del file di output per scrivere i risultati
    std::ofstream fout{"data/average.dat"};

    // Calcolo della media
    for(int i = 0; i < N; i++){
        sum = 0;
        for(int j = 0; j < L; j++){
            sum += rnd.Rannyu(); // Generazione di un numero casuale uniforme e aggiunta alla somma
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