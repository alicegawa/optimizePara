#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <mpi.h>
#include "cmaes_interface.h"

#define I_AM_ROOT_IN_MAIN (main_myid == root_process_main)
#define I_AM_ROOT_IN_SPLIT (split_myid == root_process_split)
#define I_AM_ROOT_IN_NRN (spawn_myid == root_process_spawn)

/*definition of the wrapper of MPI_Bcast (function for communication with NEURON) */
void MPI_Bcast_to_NEURON(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    /* datatype: now support MPI_DOUBLE only*/
    MPI_Bcast(&count, 1, MPI_INT, root, comm);
    MPI_Bcast(buffer, count, datatype, root, comm);
    fflush(stdout);
}/*MPI_Bcast_to_NEURON*/

int main(int argc, char **argv){
    int main_myid, main_size;
    int split_myid, split_size;
    int spawn_myid, spawn_size;
    MPI_Comm splitcomm, intercomm, parentcomm, nrn_comm;
    int root_process_main = 0, root_process_split = 0, root_process_spawn = 0;
    int color, key;
    FILE *fp;
    char specials[] = "special";
    char **neuron_argv;
    char option_mpi[] = "-mpi", option_nobanner[] = "-nobanner", HOCFILE[] = "../hocfile/main.hoc";
    int num_of_pop_per_procs;
    int num_of_procs_nrn = 8;
    double t_start, t_end;

    cmaes_t evo;
    double *arFunvals, *arFunvals_buf1, *arFunvals_buf2, *const*pop, *xfinal;
    double *initialX, *initialSigma;
    int mu = -1;
    unsigned int seed = 0;
    unsigned int dimension = 10;/* provisional definition. This must be changed based on boundary setting*/
    int max_eval = -1, max_iter = -1;
    char initfile[] = "cmaes_initials.par";
    unsigned int num_of_pop = 8;
    double *pop_sendbuf, *pop_rcvbuf;
    double *x_temp;
    int offset;
    double flg_termination = 0;
    
    int i,j;
    double util;
/* initialize MPI settings*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

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

    /*initialize for/and cmaes settings*/
    num_of_pop_per_procs = num_of_pop / num_of_procs_nrn;
    offset = num_of_pop_per_process;
    pop_sendbuf = (double *)calloc((num_of_pop + offset) * dimension, sizeof(double));
    pop_rcvbuf = (double *)malloc(num_of_pop_per_procs * dimension * sizeof(double));
    arFunvals_buf1 = (double *)calloc(num_of_pop_per_procs, sizeof(double));
    arFunvals_buf2 = (double *)calloc(num_of_pop + offset, sizeof(double));
    x_temp = (double *)malloc(dimension * sizeof(double));
    
    initialX = (double *)malloc(sizeof(double) * dimension);
    if(initialX==NULL){ printf("memory allocation error for initialX \n"); return -1;}
    initialSigma = (double *)malloc(sizeof(double) * dimension);
    if(initialSigma==NULL){ printf("memory allocation error for initialSigma \n"); return -1;}
    srand((unsigned)time(NULL));
    util = 1.0 / ((double)RAND_MAX + 2.0);
    for(i=0; i<dimension; ++i){
	initialX[i] = 3.0 + (7.0 - 3.0) * (double)(rand() + 1.0) * util;
	initialSigma[i] = (10.0 - 0.0) * 0.25;
    }
    arFunvals = cmaes_init(&evo, dimension, initialX, initialSigma, seed, num_of_pop, mu, max_eval, max_iter, initfile);
    max_iter = (int)cmaes_Get(&evo, "MaxIter");
    max_eval = (int)cmaes_Get(&evo, "maxeval");
    
    /* setting MPI_Comm_split section*/
    key = 0;
    color = main_myid;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &splitcomm);
    MPI_Comm_rank(splitcomm, &split_myid);
    MPI_Comm_size(splitcomm, &split_size);

    /*setting the argv for spawn*/
    sprintf(neuron_argv[3],"COLOR=%d",color);
    neuron_argv[5] = NULL;

    /*********caution********/
    /*probably the number of split comm is equal to the number of the main process, so the if sequence below may be unnecessary.*/
    /*********caution********/
    /*execute MPI_Comm_spawn (from below, until sending finalize signal to nrncomm, the program should not end.*/
    if( I_AM_ROOT_IN_SPLIT){
	MPI_Comm_spawn(specials, neuron_argv, 4, MPI_INFO_NULL, 0, splitcomm, &intercomm, MPI_ERRCODES_IGNORE);
	MPI_Intercomm_merge(intercomm, 0, &nrn_comm);
	MPI_Comm_size(nrn_comm, &spawn_size);
	MPI_Comm_rank(nrn_comm, &spawn_id);
    }
    /*if you want to send information to NEURON in setting, send here*/
    if( I_AM_ROOT_IN_SPLIT){
	send_count = 3;/*template*/
	double info[] = {1,0, 1.0, 1.0};/*template*/
	MPI_BCast_to_NEURON(info, send_count, MPI_DOUBLE, root_process_spawn, nrn_comm);
    }
    fflush(stdout);

    tstart = MPI_Wtime();
    while(1){
	if(I_AM_ROOT_IN_SPLIT){
	    pop = cmaes_SamplePopulation(&evo);/*do not change content of pop*/
	    for(i=0;i<num_of_pop;++i){
		/*write the boundary transformation here, under construction*/
		for(j=0;j<dimension;++j){
		    pop_sendbuf[(i + offset) * dimension + j] = x_temp[j];
		}
	    }
	    /*evaluate the new searching points*/
	    MPI_Scatter(pop_sendbuf, num_of_pop_per_procs * dimension, MPI_DOUBLE, pop_rcvbuf, num_of_pop_per_procs * dimension, MPI_DOUBLE, root_process_spawn, nrn_comm);
	    /* wait for NEURON simulation in worker nodes */
	    MPI_Gather(arFunvals_buf1, num_of_pop_per_procs, MPI_DOUBLE, arFunvals_buf2, num_of_pop_per_procs, MPI_DOUBLE, root_process_spawn, nrn_comm);

	    for(i=0;i<num_of_pop; ++i){
		arFunvals[i] = arFunvals_buf[i + offset];
	    }
	    /*update the search distribution used for cmaes_sampleDistribution()*/
	    cmaes_UpdateDistribution(&evo, arFunvals); /*assume that pop[] has not been modified*/
	}
	fflush(stdout);
	/*terminatinn*/
	if(I_AM_ROOT_IN_SPLIT){
	    if(cmaes_TestForTermination(&evo)){
		flg_termination = 1;
	    }
	    send_count = 1;
	    MPI_Bcast_to_NEURON(&flg_termination, 1, MPI_DOUBLE, root_process_split, splitcomm);
	    MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, root_process_split, splitcomm);
	    if((int)flg_termination){
		break;
	    }
	}/*cmaes termination*/
    }/*end of cmaes loop*/

    MPI_Barrier(MPI_COMM_WORLD);
    t_end = MPI_Wtime();

    if( I_AM_ROOT_IN_SPLIT){
	printf("#Stop:\n%s\n", cmaes_TestForTermination(&evo)); /*print termination reason*/
	printf("\n# operation finished.\n# elapsed time: %f\n #fbest: %f\n #xbest:", t_end - t_start, cmaes_Get(&evo, "fbestever"));
	/*exec "my boundary transformation" and print out "xbest"*/
	fflush(stdout);
    }

    free(x_temp);
    free(initialX);
    free(initialSigma);
    free(pop_sendbuf);
    free(pop_rcvbuf);
    free(arFunvals_buf1);
    free(arFunvals_buf2);
    
    for(i=0;i<6;i++){
	free(neuron_argv[i]);
    }
    free(neuron_argv);

    MPI_Finalize();

    return 0;
}
