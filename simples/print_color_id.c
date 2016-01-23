#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){
    int color;
    int numprocs, myid;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if(argc<2){
	printf("suffix is not enough\n");
	return -1;
    }

    color = atoi(argv[1]);
    printf("I am %d color\'s %d th process\n", color, myid);

    MPI_Finalize();
    return 0;
}
