#include<stdio.h>
#include<mpi.h>

int main(int argc, char **argv){
    int num;
    int numprocs, myid;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if(argc<2){
	printf("suffix is not enough\n");
	return -1;
    }

    num = atoi(argv[1]);
    printf("input number is %d\n", num);

    MPI_Finalize();
    return 0;
}
