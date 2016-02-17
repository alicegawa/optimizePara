#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){
    int i;
    int numprocs, myid, spawn_size, spawn_id;
    char command[] = "./simple";
    MPI_Comm comm1, comm2, parentcomm;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Comm_get_parent(&parentcomm);
    
    MPI_Comm_spawn(command, MPI_ARGV_NULL, 4, MPI_INFO_NULL, 0, MPI_COMM_SELF, &comm1, MPI_ERRCODES_IGNORE);
    printf("finish comm spawn\n");
    fflush(stdout);
    MPI_Comm_free(&comm1);

   
    MPI_Finalize();

    return 0;
}
