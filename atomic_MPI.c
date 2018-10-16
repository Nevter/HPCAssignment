/**
 * MPI Implementation of the atomic modeling thing
 */
#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])  {
    
    int myid;
    
    MPI_Init(&argc, &argv);                 //Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);   //get rank of node's process
    
    printf("Hello MPI World from process %d\n", myid);
    
    MPI_Finalize();

    return 0;
}
