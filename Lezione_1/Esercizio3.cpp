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

double error(double av, double av2, int n){
    return n == 0 ? 0 : sqrt((av2 - av*av)/n);
}

double length(double x,double y){
    return sqrt(x*x + y*y);
}

int main(){
    Random rnd; // Oggetto per la generazione di numeri casuali
    rnd.initialize(); // Inizializzazione del generatore di numeri casuali (seed e numeri primi)

    // Apertura del file di output per scrivere i risultati
    std::ofstream fout{"data/buffon.dat"};

    //Esperimento di Buffon

    int M{1000000}; // Numero di lanci
    int N{100}; // Numero di blocchi
    int N_thr{int(M/N)}; // Numero di lanci per blocco
    double d{1}; // Distanza tra le righe
    double L{0.5}; // Lunghezza dell'ago
    double y1{}; // Posizione y del centro dell'ago
    double x2{}; // Posizione x dell'estremita' dell'ago
    double y2{}; // Posizione y dell'estremita' dell'ago
    int N_hit{};
    double pi{};
    double sum_pi{};
    double sum_pi2{};

    for(int i = 0; i < N; i++){ // Ciclo sui blocchi
        N_hit = 0;
        for(int j = 0; j < N_thr; j++){
            y1 = rnd.Rannyu(0., d); // Genero la posizione y del centro dell'ago
            //genero punti distribuiti su una circonferenza di raggio 1
            do{ 
                x2 = rnd.Rannyu();
                y2 = rnd.Rannyu(); 
            }while(length(x2,y2) > 1); 
            y2 /= length(x2,y2); // In questo modo y2 e' distribuito come sin(theta) con theta uniformemente distribuito tra 0 e pi/2
            y2 = y1 + L*y2;

            if(y2 <=0 || y2 >= d){ // Controllo se l'ago ha colpito una riga
                N_hit++;
            }
        }
        pi = (2*L*N_thr)/(d*N_hit); //Calcolo di pi
        sum_pi += pi; // Somma parziale
        sum_pi2 += pi*pi; // Somma parziale di pi^2 
        fout << (i+1)*N_thr << " " << sum_pi/double(i+1) << " " << error(sum_pi/double(i+1), sum_pi2/double(i+1), i) << std::endl;
    }

    return 0;
}