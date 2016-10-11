#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char** argv){
    int num[4];
    int sum=0;
    int i;
    int numprocs, myid;

    usleep(10);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    
    if(argc<5){
	printf("suffix is not enough\n");
	//return -1;
	for(i=0;i<4;++i){
	    num[i] = i;
	    sum += num[i];
	}
    }else{
	for(i=0;i<4;i++){
	    num[i] = atoi(argv[i+1]);
	    sum += num[i];
	}
    }
    printf("sum = %d\n",sum);

    MPI_Finalize();
    
    return 0;
}
    
