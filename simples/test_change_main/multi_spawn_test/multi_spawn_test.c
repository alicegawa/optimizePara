#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int i;
    int main_myid, main_numproc;
    char command[] = "./call_add";
    MPI_Comm child_comm[4], parentcomm, all_spawns;
    int spawn_size[4], spawn_myid[4];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    MPI_Comm_get_parent(&parentcomm);

    if(main_myid==0){
	for(i=0; i<4; ++i){
	    MPI_Comm_spawn(command, NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_SELF, &child_comm[i], MPI_ERRCODES_IGNORE);
	    MPI_Comm_size(child_comm[i], &spawn_size[i]);
	    MPI_Comm_rank(child_comm[i], &spawn_myid[i]);
	    printf("child_comm[%d] 's size is %d\n", i, spawn_size[i]);
	}
    }

    if(main_myid == 0){
	for(i=0;i<4;++i){
	    MPI_Comm_free(&child_comm[i]);
	}
    }
    
    MPI_Finalize();
    return 0;
}
