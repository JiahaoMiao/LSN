#include "mpi.h"
#include <iostream>
using namespace std;
int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int isend[2],irecv[2];
    for(int i=0;i<2;i++)
        isend[i]=rank*2+i;
    MPI_Reduce(&isend[0],&irecv[0],1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&isend[1],&irecv[1],1,MPI_INTEGER,MPI_PROD,0,MPI_COMM_WORLD);
    if(rank==0)
        cout<<"irecv[0]="<<irecv[0]<<", irecv[1]="<<irecv[1]<<endl;

    MPI_Finalize();
    return 0;
}