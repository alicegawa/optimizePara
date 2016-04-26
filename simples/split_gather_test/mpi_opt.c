#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int n, myid, numprocs, i;
    int new_myid, new_size;
    MPI_Comm comm1, comm2, parentcomm;
    int color, key;
    FILE *fp;
    int num_memofgroup = 2, call_final;
    int maxprocs = 2;
    int *test_gather_global, *test_gather_split;
    int *test_scatter_global;
    int test_scatter_global_rcv[2];
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    test_gather_global = (int *)calloc(numprocs, sizeof(int));
    test_scatter_global = (int *)calloc(2 * numprocs, sizeof(int));
    for(i=0;i<numprocs*2;i++){
	test_scatter_global[i] = i;
    }
    
    if(argc>1){
	num_memofgroup = atoi(argv[1]);
    }
                
    key = 0;
    color = myid/num_memofgroup;
    call_final = numprocs/num_memofgroup;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm1);
    MPI_Comm_rank(comm1, &new_myid);
    MPI_Comm_size(comm1, &new_size);
    test_gather_split = (int *)calloc(new_size, sizeof(int));
    
    MPI_Gather(&new_myid, 1, MPI_INT, test_gather_split, 1, MPI_INT, 0, comm1);
    if(new_myid==0){
	printf("result of test_gather_split of %d color: %d, %d\n",color, test_gather_split[0], test_gather_split[1]);
    }
    MPI_Gather(&myid, 1, MPI_INT, test_gather_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(myid==0){
	printf("result of test_gather_global: ");
	for(i=0;i<numprocs;i++){
	    printf("%d ", test_gather_global[i]);
	}
	printf("\n");
    }
    MPI_Scatter(test_scatter_global, 2, MPI_INT, test_scatter_global_rcv, 2, MPI_INT, 0, MPI_COMM_WORLD);
    if(new_myid == 0){
	for(i=0;i<2;i++){
	    printf("result of test_scatter_global_rcv of %dth process in %d color: %d, %d\n",new_myid, color, test_scatter_global_rcv[0], test_scatter_global_rcv[1]);
	}
    }

    free(test_gather_global);
    free(test_gather_split);
    MPI_Finalize();
}
