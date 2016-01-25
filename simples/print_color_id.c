#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){
    int color;
    int numprocs, myid;
    FILE *fp;
    char filename[50];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if(argc<2){
	printf("suffix is not enough\n");
	return -1;
    }

    color = atoi(argv[1]);
    sprintf(filename,"color_%d_id_%d.txt",color, myid);
    fp = fopen(filename,"w");
    
    printf("I am %d color\'s %d th process\n", color, myid);
    fprintf(fp, "I a %d color\"s %d th process \n", color, myid);
    fclose(fp);
    MPI_Finalize();
    return 0;
}
