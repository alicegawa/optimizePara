#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv){

    MPI_Comm parent;
    int numprocs, myid;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//MPI_Comm_get_parent(&parent);
    printf("hello simple \n");
    MPI_Finalize();
    return 0;
}
