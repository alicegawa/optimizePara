#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char **argv){
    int i;
    int main_myid, main_numproc;
    int parent_size, parent_rank;
    char command[] = "./add";
    MPI_Comm child_comm[4], parentcomm, spawn_comm;
    int scatter_sendvec[5], scatter_rcvvec[5];
    int gather_sendvec[5], gather_rcvvec[6];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    MPI_Comm_get_parent(&parentcomm);

    MPI_Intercomm_merge(parentcomm, 1, &spawn_comm);
    MPI_Comm_size(spawn_comm, &parent_size);//this sentense is required. without it, program stops here...
    MPI_Comm_rank(spawn_comm, &parent_rank);
    printf("parent comm size is %d (my_rank is %d)\n", parent_size, parent_rank);

    //printf("in call_add, main_numproc is %d\n", main_numproc);

     for(i=0;i<5;++i){
	scatter_rcvvec[i] = -1;
    }

    gather_sendvec[0]= main_myid;

    MPI_Scatter(scatter_sendvec, 1, MPI_INT, scatter_rcvvec, 1, MPI_INT, 0, spawn_comm);
 
    printf("myid = %d, scatter_rcvvec[0] = %d\n", main_myid, scatter_rcvvec[0]);
    
    MPI_Gather(gather_sendvec, 1, MPI_INT, gather_rcvvec, 1, MPI_INT, 0, spawn_comm);
   
 
     MPI_Comm_spawn(command, NULL, 4, MPI_INFO_NULL, 0, MPI_COMM_SELF, &child_comm[main_myid], MPI_ERRCODES_IGNORE);

    
    MPI_Comm_free(&child_comm[main_myid]);

    MPI_Finalize();

    return 0;
}
