#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int myid, numprocs, i;
    int new_myid, new_size;
    MPI_Comm comm1, comm2, parentcomm;
    int color, key;
    int num_memofgroup=1;
    char command[]="./print_color_id";
    char **command_argv;
        
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Comm_get_parent(&parentcomm);

    command_argv = (char **)malloc(sizeof(char *) * 2);
    for(i=0;i<2;i++){
	command_argv[i] = (char *)malloc(sizeof(char) * 3);
    }
    
    if(argc>1){
	num_memofgroup = atoi(argv[1]);
    }
        
    key = 0;
    color = myid/num_memofgroup;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm1);
    MPI_Comm_rank(comm1, &new_myid);
    MPI_Comm_size(comm1, &new_size);

    sprintf(command_argv[0], "%d", color);
    command_argv[1] = NULL;
    
    MPI_Comm_spawn(command, command_argv, 8, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("finish the program\n");

    for(i=0;i<2;i++){
	free(command_argv[i]);
    }
    free(command_argv);
    
    MPI_Finalize();
}
