#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

void gather_and_scatter(int *scatter_sendvec, int *scatter_rcvvec, int scatter_sendnum, int *gather_sendvec, int*gather_rcvvec, int gather_sendnum, MPI_Comm spawn_comm){
    int i;
    int scatter_sendvec2[5];
    int gather_rcvvec2[5];
    for(i=0;i<5;++i){
	scatter_sendvec2[i] = i + 100;
    }
    MPI_Scatter(scatter_sendvec2, 1, MPI_INT, scatter_rcvvec, 1, MPI_INT, 0, spawn_comm);
    gather_sendvec[0] = -100;
    MPI_Gather(gather_sendvec, 1, MPI_INT, gather_rcvvec2, 1, MPI_INT, 0, spawn_comm);
    for(i=0;i<5;++i){
	printf("gather_rcvvec2[%d] = %d\n", i, gather_rcvvec2[i]);
    }
}

void gather_and_scatter2(int *scatter_sendvec, int *scatter_rcvvec, int scatter_sendnum, int *gather_sendvec, int*gather_rcvvec, int gather_sendnum, MPI_Comm spawn_comm){
    gather_and_scatter(scatter_sendvec, scatter_rcvvec, 1, gather_sendvec, gather_rcvvec, 1, spawn_comm);
}

void gather_and_scatter3(int scatter_sendnum, int gather_sendnum, MPI_Comm spawn_comm){
    int i;
    int *scatter_sendvec2, *gather_rcvvec2;
    int *scatter_rcvvec, *gather_sendvec;
    scatter_sendvec2 = (int *)malloc(sizeof(int) * 5);
    gather_rcvvec2 = (int *)malloc(sizeof(int) * 5);
    scatter_rcvvec = (int *)malloc(sizeof(int));
    gather_sendvec = (int *)malloc(sizeof(int));
    for(i=0;i<5;++i){
	scatter_sendvec2[i] = i + 100;
    }
    MPI_Scatter(scatter_sendvec2, 1, MPI_INT, scatter_rcvvec, 1, MPI_INT, 0, spawn_comm);
    gather_sendvec[0] = 100;
    MPI_Gather(gather_sendvec, 1, MPI_INT, gather_rcvvec2, 1, MPI_INT, 0, spawn_comm);
    for(i=0;i<5;++i){
	printf("gather_rcvvec2[%d] = %d\n", i, gather_rcvvec2[i]);
    }
    free(scatter_sendvec2); scatter_sendvec2 = NULL;
    free(gather_rcvvec2); gather_rcvvec2 = NULL;
}

int main(int argc, char **argv){
    int i;
    int main_myid, main_numproc;
    char command[] = "./call_add";
    char *command2[] = {"./call_add", "./call_add", "./call_add", "./call_add"};
    int array_of_maxprocs[] = {1, 1, 1, 1};
    MPI_Info array_of_info[] = {MPI_INFO_NULL, MPI_INFO_NULL, MPI_INFO_NULL, MPI_INFO_NULL};
    MPI_Comm spawn_comm, parentcomm, intercomm;
    int spawn_size, spawn_myid;
    int scatter_sendvec[5] = {0, 1, 2, 3, 4};
    int scatter_rcvvec[5];
    int gather_sendvec[5], gather_rcvvec[6];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    MPI_Comm_get_parent(&parentcomm);

    /* if(main_myid==0){ */
    /* 	for(i=0; i<4; ++i){ */
    /* 	    MPI_Comm_spawn(command, NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_SELF, &child_comm[i], MPI_ERRCODES_IGNORE); */
    /* 	} */
    /* } */

    /* MPI_Comm_spawn_multiple(4, command2, NULL, array_of_maxprocs, array_of_info, 0, MPI_COMM_SELF, &child_comm[0], MPI_ERRCODES_IGNORE); */

    MPI_Comm_spawn(command, NULL, 4, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(intercomm, 0, &spawn_comm);
    MPI_Comm_size(spawn_comm, &spawn_size);
    MPI_Comm_rank(spawn_comm, &spawn_myid);
    printf("spawn_comm 's size is %d\n", spawn_size);

    /* MPI_Scatter(scatter_sendvec, 1, MPI_INT, scatter_rcvvec, 1, MPI_INT, 0, spawn_comm); */

    /* gather_sendvec[0] = -100; */

    /* MPI_Gather(gather_sendvec, 1, MPI_INT, gather_rcvvec, 1, MPI_INT, 0, spawn_comm); */

    //gather_and_scatter2(scatter_sendvec, scatter_rcvvec, 1, gather_sendvec, gather_rcvvec, 1, spawn_comm);
    gather_and_scatter3(1, 1, spawn_comm);

    for(i=0;i<5;++i){
    	//printf("gather_rcvvec[%d] = %d\n", i, gather_rcvvec[i]);
    }

    /* if(main_myid == 0){ */
    /* 	for(i=0;i<4;++i){ */
    /* 	    MPI_Comm_free(&child_comm[i]); */
    /* 	} */
    /* } */
    MPI_Comm_free(&intercomm);
    MPI_Comm_free(&spawn_comm);
    

    MPI_Finalize();
    return 0;
}
