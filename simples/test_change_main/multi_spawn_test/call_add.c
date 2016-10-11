#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char **argv){
    int i;
    int main_myid, main_numproc;
    char command[] = "./add";
    MPI_Comm child_comm[4], parentcomm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    MPI_Comm_get_parent(&parentcomm);

    printf("in call_add, main_numproc is %d\n", main_numproc);

    MPI_Comm_spawn(command, NULL, 4, MPI_INFO_NULL, 0, MPI_COMM_SELF, &child_comm[main_myid], MPI_ERRCODES_IGNORE);
    
    MPI_Comm_free(&child_comm[main_myid]);

    MPI_Finalize();

    return 0;
}
