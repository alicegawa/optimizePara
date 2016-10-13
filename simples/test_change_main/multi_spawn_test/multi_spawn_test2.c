#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

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

    MPI_Scatter(scatter_sendvec, 1, MPI_INT, scatter_rcvvec, 1, MPI_INT, 0, spawn_comm);

    gather_sendvec[0] = -100;

    MPI_Gather(gather_sendvec, 1, MPI_INT, gather_rcvvec, 1, MPI_INT, 0, spawn_comm);

    for(i=0;i<5;++i){
    	printf("gather_rcvvec[%d] = %d\n", i, gather_rcvvec[i]);
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
