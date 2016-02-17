#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int myid, numprocs, i;
    int split_myid, split_size;
    MPI_Comm splitcomm, intercomm, parentcomm, nrn_comm;
    int color, key;
    FILE *fp;
    char specials[] = "special";
    char **neuron_argv;
    char option_mpi[] = "-mpi", option_nobanner[] = "-nobanner", HOCFILE[] = "../hocfile/main.hoc";
    

    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Comm_get_parent(&parentcomm);

    neuron_argv = (char **)malloc(sizeof(char *) * 6);
    for(i=0;i<6;i++){
	neuron_argv[i] = (char *)malloc(sizeof(char) * 10);
	if(i==0){
	    sprintf(neuron_argv[i],"%s",option_mpi);
	}else if(i==1){
	    sprintf(neuron_argv[i],"%s",option_nobanner);
	}else if(i==2){
	    sprintf(neuron_argv[i],"-c");
	}else if(i==4){
	    sprintf(neuron_argv[i],"%s",HOCFILE);
	}else{
	}
    }

    key = 0;
    color = myid;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &splitcomm);
    MPI_Comm_rank(splitcomm, &split_myid);
    MPI_Comm_size(splitcomm, &split_size);

    sprintf(neuron_argv[3],"COLOR=%d",color);
    neuron_argv[5] = NULL;

    MPI_Comm_spawn(specials, neuron_argv, 4, MPI_INFO_NULL, 0, splitcomm, &intercomm, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(intercomm, 0, &nrn_comm);

    fflush(stdout);

    for(i=0;i<6;i++){
	free(neuron_argv[i]);
    }
    free(neuron_argv);

    MPI_Finalize();

    return 0;
}
