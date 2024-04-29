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
    std::ofstream fout{"data/Chi2.dat"};

    // Parametri per il calcolo del chi quadro
    int M{100}; // Numero di bins 
    int n{10000}; // Numero di estrazioni casuali 
    double expected = double(n) / double(M); // Numero di eventi attesi per ogni bin
    int chi[M]; // Array per memorizzare il numero di eventi osservati per ogni bin
    for(int i = 0; i < M; i++) chi[i] = 0; // Inizializzazione dell'array chi

    double x{}, chi2{};
    int bin_index{};

    for(int j=0; j < 10000; j++){ // Ciclo 10000 volte
        // Generazione di numeri casuali e calcolo del chi quadro
        for(int i = 0; i < n; i++){
            x = rnd.Rannyu(); // Generazione di un numero casuale
            bin_index = int(x * M); // Calcolo del bin in cui cade il numero casuale
            chi[bin_index]++; // Incremento del numero di eventi osservati per il bin corrispondente
        }

        // Calcolo del chi quadro
        chi2  = 0;
        for(int i = 0; i < M; i++){
            chi2 += std::pow(chi[i] - expected, 2) / expected; // Calcolo del chi quadro
            chi[i] = 0; // Reset del numero di eventi osservati per ogni bin
        }
        fout << chi2 << std::endl; // Scrittura del risultato su file
    }

    return 0;
}   