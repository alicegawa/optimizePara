#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int n, myid, numprocs, i;
    int new_myid, new_size;
    MPI_Comm comm1, comm2, parentcomm, nrn_comm;
    int color, key;
    FILE *fp;
    char filename[100];
    int data[4];
    int num_memofgroup=1, call_final;
    char command[]="./print_color_id";
    char specials[]="./special";
    char *neuron_argv[] = {"-mpi", "-nobanner", "para_test.hoc", NULL};
    char **command_argv;
    char option_mpi[] ="-mpi", HOCFILE[] = "para_test.hoc", option_nobanner[]="-nobanner";
    int maxprocs=2;
    
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
        
    for(i=0;i<4;i++){
	data[i] = 0;
    }
        
    key = 0;
    color = myid/num_memofgroup;
    call_final = numprocs/num_memofgroup;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm1);
    MPI_Comm_rank(comm1, &new_myid);
    MPI_Comm_size(comm1, &new_size);

    sprintf(command_argv[0], "%d", color);
    command_argv[1] = NULL;
    
    MPI_Comm_spawn(command, command_argv, 4, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    fflush(stdout);

    for(i=0;i<2;i++){
	free(command_argv[i]);
    }
    free(command_argv);
    
    MPI_Finalize();
}
