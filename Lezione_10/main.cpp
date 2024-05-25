#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "RndGen/random.h"
#include "TSP.h"

using namespace std;

void showProgressBar(int current, int total) {
    int barWidth = 70;
    float progress = (float)current / (float)total;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}


int main(int argc, char* argv[]){

    //parameters
    int L = 3000;//number of paths
    int M = 10000;//number of generations
    bool migration = true;
    int Mmigr = 30;//number of generations between migrations
    

    int size{},rank{};

    // Start MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);//how many of us?
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);//who am i? 
    MPI_Status stat; //status of the communication

    //each process will have its own random generator
    Random rnd;
    rnd.initialize(rank);

    std::vector<city> cities;
    string filename{"DATA/"};

    if(migration){
        filename += "Migration/";
    }else{
        filename += "NoMigration/";
    }
    filename += to_string(rank)+"/";

    //Common for all processes (maybe not the best way to do it but it works for now)
    ifstream Input("INPUT/cap_prov_ita.dat");
    if(!Input){
        cerr << "Error opening the input file" << "\n";
        return -1;
    }
    int N{};
    double x{},y{};
    while(Input >> x >>y){
        cities.emplace_back(city(x,y));
        N++;
    }
    cities.reserve(N);
    cities.shrink_to_fit();
    Input.close();

    //Print initial cities
    if(rank ==0){ //Only the master will do this
        ofstream config(filename+"config.dat");
        for(auto c : cities){
            config << c.getX() << " " << c.getY() << "\n";
        }
        //to close the path
        config << cities[0].getX() << " " << cities[0].getY() << "\n";
        config.close();    
    }
    
    //start with the ordered path
    //path always starts from 0 so i omit it and work with the rest of the cities
    vector<int> way(N-1);
    for(int i = 1; i < N; i++){
        way[i-1] = i;
    }
    path route{way};
    vector<path> paths(L);
    int i{},j{};
    //generate L paths with random swaps 
    for(int k{};k<L;k++){
        for(int l{};l<3*N;l++){//2*N swaps
            i = rnd.Rannyu(0,N-1);
            j = rnd.Rannyu(0,N-1);
            while(i==j){
                j= rnd.Rannyu(0,N-1);
            }
            // i % (N-1) and j % (N-1) to avoid the value N-1
            route.swap(i%(N-1),j%(N-1));//swap of two cities
        }
        paths[k] = route;
    }
    TSP tsp{cities,paths,rnd};
     
    // std::cout << "Process " << rank << " Best Initial path: " << "\n";
    // tsp.PrintBest();

    ofstream ave{filename+"AverageLength.dat"};
    ofstream best{filename+"BestLength.dat"};
    if(!ave || !best){
        cerr << "Error opening the output file" << "\n";
        return -1;
    }

    ave << 0 << " " << tsp.AverageLength() << "\n";

    //for the gif, not needed now
    //tsp.BestPath(filename+"config/BestPath0.dat");

    for(int i{1};i<=M;i++){
        if(!migration){
            tsp.NewGeneration();  //no migration
        }else{//with migration
            if(i%Mmigr==0){
                //first migrate then generate new generation
                int giver{};
                int receiver{};
                vector<int> migrator(N-1); //path without the first city
                int itag=1;

                // only the master will choose the individual to migrate
                if (rank == 0) {
                    giver = static_cast<int>(rnd.Rannyu(0, size)); // choose the starting core
                    receiver = static_cast<int>(rnd.Rannyu(0, size)); // choose the destination core
                    while (receiver == giver) {
                        receiver = static_cast<int>(rnd.Rannyu(0, size)); // ensure receiver is different from giver
                    }
                    std::cout << "\nMigration gen: " << i << " (giver, receiver): (" << giver << ", " << receiver << ")\n";
                }
                // master broadcasts the migration data
                MPI_Bcast(&giver, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&receiver, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // if I am the giver I send the best individual
                if(rank == giver){
                    for(int i = 0; i < N - 1; i++){
                        migrator[i] = tsp[0][i];
                    }
                    MPI_Send(migrator.data(), N - 1, MPI_INT, receiver, itag, MPI_COMM_WORLD);
                }
                // if I am the receiver I receive the individual and substitute the worst individual
                if(rank == receiver){
                    MPI_Recv(migrator.data(), N - 1, MPI_INT, giver, itag, MPI_COMM_WORLD, &stat);
                    tsp[L-1] = path(migrator); // move assignment
                }
            }

            tsp.NewGeneration();
        }

        ave << i << " " << tsp.AverageLength() << "\n";
        best << i << " " << tsp.BestLength() << "\n";
        // if(i%2==0){
        // tsp.BestPath(filename+"config/BestPath"+to_string(i)+".dat");
        // }
        // Show progress bar
        if(rank == 0){
            showProgressBar(i, M);
        }
    }
    //print the best path to file
    tsp.BestPath(filename+"BestPath.dat");
    MPI_Finalize();
    
    return 0;
}