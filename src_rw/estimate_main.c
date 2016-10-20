#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "cmaes_interface.h"
#include "my_boundary_transformation.h"

#define I_AM_ROOT_IN_MAIN (main_myid == root_process_main)
#define I_AM_ROOT_IN_SPAWN (split_myid == root_process_spawn)

#define TEXT_BUFFER_SIZE 2048
#define MAX_NUM_PARAM 150
#define MAX_NUM_TARGET 5

/* under construction*/
int loadRangeFile(char *filename, my_boundary_transformation_t *t){
    int dimension = 0;
    char buf[TEXT_BUFFER_SIZE];
    double lowerBounds[MAX_NUM_PARAM];
    double upperBounds[MAX_NUM_PARAM];
    unsigned int flg_log[MAX_NUM_PARAM];

    FILE *fp;
    if((fp=fopen(filename, "r"))==NULL){
	printf("file open eror\n");
	exit(EXIT_FAILURE);
    }
    while( fgets(buf, TEXT_BUFFER_SIZE, fp) != NULL){
	if(strncmp(buf, "#", 1) == 0){ continue;}
	sscanf(buf, "%*s\t%lf\t%lf\t%*lf\t%d\n", &lowerBounds[dimension], &upperBounds[dimension], &flg_log[dimension]);
	//flg_log[dimension] = (unsigned char)atoi((const char*)&flg_log[dimension]);
	dimension++;//make the dimension information here
    }
    my_boundary_transformation_init(t, lowerBounds, upperBounds, flg_log, dimension);
    fclose(fp);
    return 0;
}

/*print gene data to file*/
int printGene(FILE *fp, const double *x, int dimension){
    int i;
    for(i=0; i<dimension; ++i){
	fprintf(fp, "%f\t", x[i]);
    }
    return 0;
}/*printGene()*/

int main(int argc, char **argv){
    int main_myid, main_size;
    int spawn_myid, spawn_size;
    MPI_Comm parentcomm, intercomm, spawn_comm;
    int root_process_main = 0, root_process_spawn = 0;
    
    FILE *fp;
    char make_neuro_spawn[] = "./make_neuro_spawn";
    char **spawn_argv;//information that is passed to spawn procs
    int spawn_argv_size = 3;/*temporary setting*/

    unsigned int num_of_pop = 2, num_of_pop_per_child, num_of_child_procs = 2;
    cmaes_t evo;
    my_boundary_transformation_t my_boundaries;
    double *arFunvals, *arFunvals_whole, *arFunvals_whole_buf;
    double *const*pop;
    double *xfinal, *initialX, *initialSigma, *x_temp, *x_temp_temp;
    int mu = -1;
    unsigned int seed = 0, dimension;
    int max_iter = -1, max_eval = -1;
    char initfile[] = "./cmaes_initials.par";
    
    double *pop_sendbuf_spawn_whole, *pop_sendbuf_spawn_weight, *pop_sendbuf_spawn_delay;
    double *pop_rcvbuf_spawn_whole, *pop_rcvbuf_spawn_weight, *pop_rcvbuf_spawn_delay;
    int offset_gather, offset_scatter;
    int num_of_weights, num_of_delays, num_of_sendparams_per_child, num_of_sendparams_per_child_sep;
    double flg_termination = 0;
    char range_filename[] = "../data/params.txt";
    
    int i, j, k;
    double util;
    int loop_count = 0;

    int n_run = 0;
    int countevals, generation;
    int gen_restart[] = { 1, 2, 4, 8, 16, 32, 64};
    double *restartX, *restartSigma;
    double restartSigma_defaults = 2.0;

    double t_start, t_end;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);
    
    MPI_Comm_get_parent(&parentcomm);

    /* check input params */
    if(argc >= 2){num_of_pop = atoi(argv[1]);}
    if(argc >= 3){max_iter = atoi(argv[2]);}
    if(argc >= 4){max_eval = atoi(argv[3]);}
    if(argc >= 5){num_of_child_procs = atoi(argv[4]);}
    if(argc >= 6){mu = atoi(argv[5]);}

    if(num_of_pop%num_of_child_procs != 0){
	num_of_pop -= num_of_pop%num_of_child_procs;
    }

    /* allocate spawn_argv array */
    spawn_argv = (char **)malloc(sizeof(char *) * spawn_argv_size);
    
    for(i=0;i<spawn_argv_size;++i){
	spawn_argv[i] = (char *)malloc(sizeof(char) * 256);
    }/*add setting later*/
    spawn_argv[spawn_argv_size-1] = NULL;
    
    /* configuration of the max num of iteration and evaluation */
    if(max_iter == -1){
	if(max_eval == -1){
	    max_iter = 2;
	    max_eval = max_iter * num_of_pop;
	}else{
	    max_iter = ceil((double)max_eval / num_of_pop);
	}
    }else{
	if(max_eval == -1){
	    max_eval = max_iter * num_of_pop;
	}else{
	    if(max_iter * num_of_pop < max_eval){
		max_eval = max_iter * num_of_pop;
	    }else{
		max_iter = ceil((double)max_eval / num_of_pop);
	    }
	}
    }
    /* setting of mu */
    if(mu==-1){
	mu = num_of_pop / 2;
    }

    /* initialize for/and cmaes settings and spawn parameters*/
    loadRangeFile(range_filename, &my_boundaries);
    dimension = my_boundaries.dimension;
    num_of_weights =  num_of_delays = dimension * num_of_pop / 2;
    num_of_pop_per_child = num_of_pop / num_of_child_procs;
    num_of_sendparams_per_child = num_of_pop_per_child * dimension; 
    num_of_sendparams_per_child_sep = num_of_sendparams_per_child / 2;
    offset_gather = num_of_pop_per_child;
    offset_scatter = num_of_sendparams_per_child; //avoid confusion, it had better that to delete this setting
    
    /*use only (1) or (2) */
    /* (1) */
    pop_sendbuf_spawn_whole = (double *)calloc(dimension * num_of_pop + num_of_sendparams_per_child, sizeof(double));
    pop_rcvbuf_spawn_whole = (double *)malloc(sizeof(double) * num_of_sendparams_per_child);
    /* (2) */
    pop_sendbuf_spawn_weight = (double *)calloc(dimension * num_of_pop / 2 + num_of_sendparams_per_child_sep, sizeof(double));
    pop_sendbuf_spawn_delay = (double *)calloc(dimension * num_of_pop / 2 + num_of_sendparams_per_child_sep, sizeof(double));
    pop_rcvbuf_spawn_weight = (double *)malloc(sizeof(double) * num_of_sendparams_per_child_sep);
    pop_rcvbuf_spawn_delay = (double *)malloc(sizeof(double) * num_of_sendparams_per_child_sep);
    
    arFunvals_whole_buf = (double *)calloc(num_of_pop_per_child, sizeof(double));
    arFunvals_whole = (double *)calloc(num_of_pop + num_of_pop_per_child, sizeof(double));
    x_temp = (double *)calloc(dimension, sizeof(double));
    initialX = (double *)malloc(sizeof(double) * dimension);
    initialSigma = (double *)malloc(sizeof(double) * dimension);
    restartSigma = (double *)malloc(sizeof(double) * dimension);
    srand((unsigned)time(NULL));
    util = 1.0 / ((double)RAND_MAX + 2.0);
   
    for(i=0; i<dimension; ++i){
	initialX[i] = 3.0 + (7.0 - 3.0) * (double)(rand() + 1.0) * util;
	initialSigma[i] = (10.0 - 0.0) * 0.25;
	restartSigma[i] = restartSigma_defaults;
    }

    /* start CMA-ES algorithm */
    arFunvals = cmaes_init(&evo, dimension, initialX, initialSigma, seed, num_of_pop, mu, max_eval, max_iter, initfile);
    max_iter = (int)cmaes_Get(&evo, "MaxIter");
    max_eval = (int)cmaes_Get(&evo, "maxeval");

    /* setting the spawn_argv */
    sprintf(spawn_argv[0], "%d", num_of_pop_per_child);
    sprintf(spawn_argv[1], "%d", dimension);

    /* make spawns */
    MPI_Comm_spawn(make_neuro_spawn, spawn_argv, num_of_child_procs, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(intercomm, 0, &spawn_comm);
    MPI_Comm_size(spawn_comm, &spawn_size);
    MPI_Comm_rank(spawn_comm, &spawn_myid);

    t_start = MPI_Wtime();

    /* start estimation*/
    while(1){
	pop = cmaes_SamplePopulation(&evo);/*do not change content of pop*/
	for(i=0;i<num_of_pop;++i){
	    my_boundary_transformation(&my_boundaries, pop[i], x_temp, main_myid);
	    for(j=0;j<dimension;++j){
		pop_sendbuf_spawn_whole[i * dimension + j] = x_temp[j];
	    }
	}
	for(i=0; i<num_of_pop; ++i){
	    for(j=0; j<(dimension/2); ++j){
		pop_sendbuf_spawn_weight[i * dimension / 2 + j + num_of_sendparams_per_child_sep] = pop_sendbuf_spawn_whole[i * dimension + j];
		pop_sendbuf_spawn_delay[i * dimension / 2 + j + num_of_sendparams_per_child_sep] = pop_sendbuf_spawn_whole[i * dimension + j + dimension / 2];
	    }
	}

	/* use only (1) or (2) */
	/*(1)*/
	//MPI_Scatter(pop_sendbuf_spawn_whole, num_of_sendparams_per_child, MPI_DOUBLE, pop_rcvbuf_spawn_whole, num_of_sendparams_per_child, MPI_DOUBLE, root_process_spawn, spawn_comm);
	/*(2)*/
	MPI_Scatter(pop_sendbuf_spawn_weight, num_of_sendparams_per_child_sep, MPI_DOUBLE, pop_rcvbuf_spawn_weight, num_of_sendparams_per_child_sep, MPI_DOUBLE, root_process_spawn, spawn_comm);
	MPI_Scatter(pop_sendbuf_spawn_delay, num_of_sendparams_per_child_sep, MPI_DOUBLE, pop_rcvbuf_spawn_delay, num_of_sendparams_per_child_sep, MPI_DOUBLE, root_process_spawn, spawn_comm);
    
	/* wait for calculation of child and grandchild processes */
	MPI_Gather(arFunvals_whole_buf, num_of_pop_per_child, MPI_DOUBLE, arFunvals_whole, num_of_pop_per_child, MPI_DOUBLE, root_process_spawn, spawn_comm);

	/* register the arFunvals_whole to CMA-ES's arFunval */
	for(i=0;i<num_of_pop;++i){
	    arFunvals[i] = arFunvals_whole[i + num_of_pop_per_child];
	    if(arFunvals[i] == 0){
		printf("it may be communication error\n");
		//break;
	    }
	}
	/* update CMA-ES */
	cmaes_UpdateDistribution(&evo, arFunvals);

	/* check the condition of termination */
	if(cmaes_TestForTermination(&evo)){
	    flg_termination = 1;
	}
	MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, root_process_spawn, spawn_comm);
	if((int)flg_termination){
	    break;
	}
	
	/* restart */
	if((int)cmaes_Get(&evo, "generation") == gen_restart[n_run] && I_AM_ROOT_IN_MAIN){
	    generation = cmaes_Get(&evo, "generation");
	    countevals = cmaes_Get(&evo, "eval");
	    restartX = (double *)cmaes_GetPtr(&evo, "xmean");
	    x_temp_temp = (double *)cmaes_GetPtr(&evo, "xbestever");

	    for(i=0;i<=dimension;++i){
		x_temp[i] = x_temp_temp[i];
	    }
	    arFunvals = cmaes_init(&evo, dimension, restartX, restartSigma, seed, num_of_pop, mu, max_eval, max_iter, initfile);
	    
	    evo.countevals = countevals;
	    evo.gen = generation + 1;
	 
	    for(i=0; i<=dimension; ++i){
		evo.rgxbestever[i] = x_temp[i];
		evo.rgrgx[evo.index[0]][i] = x_temp[i];
	    }
	    n_run++;
	}
	
	++loop_count;
	printf("%d times loop of cmaes finished\n", loop_count);

    }/*end of cmaes loop*/

    t_end = MPI_Wtime();

    printf("#Stop:\n%s\n", cmaes_TestForTermination(&evo)); /*print termination reason*/
    printf("\n# operation finished.\n# elapsed time: %f\n #fbest: %f\n #xbest:", t_end - t_start, cmaes_Get(&evo, "fbestever"));
    my_boundary_transformation(&my_boundaries, (double *)cmaes_GetPtr(&evo, "xbestever"), x_temp, main_myid);
    printGene(stdout, x_temp, dimension);
    printf("\n");
    fflush(stdout);

    cmaes_exit(&evo);
    my_boundary_transformation_exit(&my_boundaries);
    free(pop_sendbuf_spawn_whole);
    free(pop_rcvbuf_spawn_whole);
    free(pop_sendbuf_spawn_weight);
    free(pop_sendbuf_spawn_delay);
    free(pop_rcvbuf_spawn_weight);
    free(pop_rcvbuf_spawn_delay);
    free(arFunvals_whole_buf);
    free(arFunvals_whole);
    free(x_temp);
    free(initialX);
    free(initialSigma);
    free(restartSigma);
    
    for(i=0;i<spawn_argv_size;++i){
	free(spawn_argv[i]);
    }
    free(spawn_argv);

    MPI_Finalize();
    return 0;
}
