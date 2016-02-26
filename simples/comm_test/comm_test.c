#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


void MPI_Bcast_to_NEURON(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    MPI_Bcast(&count, 1, MPI_INT, root, comm);
    MPI_Bcast(buffer, count, datatype, root, comm);
    fflush(stdout);
}

int main(int argc, char** argv){
    int main_myid, main_size;
    int split_myid, split_size;
    int spawn_myid, spawn_size;

    MPI_Comm splitcomm, intercomm, parentcomm, nrn_comm;
    int root_process_main = 0, root_process_split = 0, root_process_spawn = 0;
    int color, key;

    FILE *fp;
    char specials[] = "special";
    char *option_mpi[] = {"-mpi","-nobanner","./para_test.hoc",NULL};
    
    int i, j;

    double *scatter1, *scatter2;
    double *gather1, *gather2;
    int send_count;
    
    double flg_termination;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    MPI_Comm_get_parent(&parentcomm);

    key = 0;
    color = main_myid;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &splitcomm);
    MPI_Comm_rank(splitcomm, &split_myid);
    MPI_Comm_size(splitcomm, &split_size);

    if(split_myid == 0){
	MPI_Comm_spawn(specials, option_mpi, 4, MPI_INFO_NULL, 0, splitcomm, &intercomm, MPI_ERRCODES_IGNORE);
	MPI_Intercomm_merge(intercomm, 0, &nrn_comm);
	MPI_Comm_size(nrn_comm, &spawn_size);
	MPI_Comm_rank(nrn_comm, &spawn_myid);
    }

   
    scatter1 = (double *)calloc(20, sizeof(double));
    scatter2 = (double *)malloc(sizeof(double)*4);
    gather1 = (double *)calloc(4, sizeof(double));
    gather2 = (double *)calloc(20, sizeof(double));
    
    while(1){
	if(split_myid==0){
	    for(i=0;i<20;i++){
		scatter1[i] = i + 10;
	    }
	    MPI_Scatter(scatter1, 4, MPI_DOUBLE, scatter2, 4, MPI_DOUBLE, 0, nrn_comm);
	    printf("end of mpi scatter\n");
	    MPI_Gather(gather1, 4, MPI_DOUBLE, gather2, 4, MPI_DOUBLE, 0, nrn_comm);
	    for(i=0;i<20;i++){
		printf("gather2[%d]=%lf\n",i, gather2[i]);
	    }
	    printf("end of mpi gather\n");
	}
	if(split_myid==0){
	    flg_termination = 1;
	    send_count = 1;
	    MPI_Bcast_to_NEURON(&flg_termination, 1, MPI_DOUBLE, root_process_spawn, nrn_comm);
	}
	MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, root_process_main, splitcomm);
	if((int)flg_termination){
	    break;
	}
    }
    printf("end of loop\n");
    fflush(stdout);
    
    MPI_Barrier(MPI_COMM_WORLD);

    free(scatter1);
    free(scatter2);
    free(gather1);
    free(gather2);

    MPI_Finalize();
    return 0;
}
