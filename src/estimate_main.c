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
#define I_AM_ROOT_IN_SPLIT (split_myid == root_process_split)
#define I_AM_ROOT_IN_NRN (spawn_myid == root_process_spawn)

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

/*definition of the wrapper of MPI_Bcast (function for communication with NEURON) */
void MPI_Bcast_to_NEURON(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    /* datatype: now support MPI_DOUBLE only*/
    MPI_Bcast(&count, 1, MPI_INT, root, comm);
    MPI_Bcast(buffer, count, datatype, root, comm);
    fflush(stdout);
}/*MPI_Bcast_to_NEURON*/

/*print gene data to file*/
int printGene(FILE *fp, const double *x, int dimension){
    int i;
    for(i=0; i<dimension; ++i){
	fprintf(fp, "%f\t", x[i]);
    }
    return 0;
}/*printGene()*/

static void Sort_score_index(double *rgFunVal, int *iindex, int n){
      int i, j;
  for (i=1, iindex[0]=0; i<n; ++i) {
    for (j=i; j>0; --j) {
      if (rgFunVal[iindex[j-1]] < rgFunVal[i])
        break;
      iindex[j] = iindex[j-1]; /* shift up */
    }
    iindex[j] = i; /* insert i */
  }
}/*Sort_score_index*/


int main(int argc, char **argv){
    int main_myid, main_size;
    int split_myid, split_size;
    int spawn_myid, spawn_size;
    MPI_Comm splitcomm, intercomm, parentcomm, nrn_comm, firstTimeWorld;
    int root_process_main = 0, root_process_split = 0, root_process_spawn = 0;
    int color, key;
    FILE *fp;
    char specials[] = "special";
    char **neuron_argv;
    char option_mpi[] = "-mpi", option_nobanner[] = "-nobanner", HOCFILE[] = "../hocfile_detail/main.hoc";
    //char *neuron_argv[] = {"-mpi", "-nobanner", "-c", "{}", "-c", "{}", "../hocfile/main.hoc", NULL};
    int num_of_pop_per_procs, num_of_pop_per_split;
    int num_of_procs_nrn = 8;//for test
    double t_start, t_end;
    int send_count;
    
    cmaes_t evo;
    double *arFunvals, *arFunvals_split_buf1, *arFunvals_split_buf2, *const*pop, *xfinal;
    double *arFunvals_whole, *arFunvals_whole_buf;
    double *initialX, *initialSigma;
    int mu = -1;//64;//-1;
    unsigned int seed = 0;
    unsigned int dimension;
    unsigned int dimension_per_nrnprocs;
    int max_eval = -1, max_iter = -1;
    char initfile[] = "cmaes_initials.par";
    unsigned int num_of_pop = 4;//1024;//for test
    double *pop_sendbuf_nrn_weight, *pop_sendbuf_nrn_delay, *pop_rcvbuf_nrn_weight, *pop_rcvbuf_nrn_delay;
    double *pop_sendbuf_split_whole, *pop_sendbuf_split_weight, *pop_sendbuf_split_delay, *pop_rcvbuf_split_weight, *pop_rcvbuf_split_delay;

    double *pop_split_whole;
    
    int num_sendparams;
    double *x_temp;
    double *x_temp_temp;
    int offset, offset_split_scatter;
    int num_params_only_weight_or_delay;
    double flg_termination = 0;
    my_boundary_transformation_t my_boundaries;
    char range_filename[] = "../data/params.txt";
    int num_of_params_per_nrnprocs;
    int num_of_weight_delay_per_procs;
    int num_of_one_gene_weight_or_delay;
    int neuron_argv_size = 6, network_size = 0;//please set in the command line
    
    int i,j,k;
    double util;
    
    int loop_count=0;
    
    /* for restart strategy*/
    int n_run=0;
    int countevals, generation;
    int gen_restart[] = { 1, 2, 3, 4, 50, 100, 150};/* temporary setting */
    double *restartX, *restartSigma;
    double restartSigma_defaults = 2.0;

    /* distributed pop update*/
    double *rgD, *x_mean, sigma, sigmasquare_div; 
    random_t rand_box;
    double const *pop_dist;
    double *sp_weight;
    int *arFunval_rank;
    double *sum_eachprocs, *sum_for_cov;
    double *sum_reduce, *sum_for_cov_reduce;/*only use in root_process_main*/
    double divider1, divider2;
    int dim_cov;
    int sum_cov_counter=0;
    
    /*test variables*/
//    double *test_sendbuf, *test_rcvbuf;
//    double *test_arFunval1, *test_arFunval2;
    /* initialize MPI settings*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);
    
    MPI_Comm_get_parent(&parentcomm);

    /*check input params */
    if(argc >= 2){num_of_pop = atoi(argv[1]);}
    if(argc >= 3){max_iter = atoi(argv[2]);}
    if(argc >= 4){max_eval = atoi(argv[3]);}
    if(argc >= 5){num_of_procs_nrn = atoi(argv[4]);}
    if(argc >= 6){mu = atoi(argv[5]);}
    if(argc >= 7){neuron_argv_size = atoi(argv[6]);}
    if(argc >= 8){network_size = atoi(argv[7]);}

    neuron_argv_size = 8;
    
    neuron_argv = (char **)malloc(sizeof(char *) * neuron_argv_size);

    
    
    for(i=0;i<neuron_argv_size;i++){
    	neuron_argv[i] = (char *)malloc(sizeof(char) * 256);
    	if(i==0){
    	    sprintf(neuron_argv[i],"%s",option_mpi);
    	}else if(i==1){
    	    sprintf(neuron_argv[i],"%s",option_nobanner);
    	}else if(i==2 || i==4){
    	    sprintf(neuron_argv[i],"-c");
    	}else if(i==(neuron_argv_size-2)){
    	    sprintf(neuron_argv[i],"%s",HOCFILE);
    	    //sprintf(neuron_argv[i],"../hocfile/main.hoc");
    	}else{
    	}
    }

    
    if(max_iter == -1){
	if(max_eval == -1){
	    max_iter = 1;
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

    if(mu==-1){
	mu = num_of_pop / 2;
    }
    
    /*initialize for/and cmaes settings*/

    loadRangeFile(range_filename, &my_boundaries);//you must modify the filename of setting file
    dimension = my_boundaries.dimension;

    
    dimension_per_nrnprocs = dimension / num_of_procs_nrn; // param dimension of each nrn process(weight and delay)
    num_params_only_weight_or_delay = dimension * num_of_pop / 2;
    
    num_of_pop_per_split = num_of_pop / main_size;
    num_of_params_per_nrnprocs = dimension_per_nrnprocs * num_of_pop_per_split;
    num_of_weight_delay_per_procs = num_of_params_per_nrnprocs / 2;

    num_sendparams = dimension * num_of_pop_per_split / 2;

    num_of_one_gene_weight_or_delay = dimension / (num_of_procs_nrn * 2);
//num_of_pop_per_procs = num_of_pop_per_split / num_of_procs_nrn;
    //offset = num_of_pop_per_procs;
    offset_split_scatter = num_of_pop_per_split;
    offset = num_of_params_per_nrnprocs / 2;

    if(I_AM_ROOT_IN_MAIN){
	pop_sendbuf_split_whole = (double *)calloc(dimension * num_of_pop, sizeof(double));
	//pop_sendbuf_split_weight = (double *)calloc(dimension * (num_of_pop + num_of_pop_per_split) / 2, sizeof(double));
	//pop_sendbuf_split_delay = (double *)calloc(dimension * (num_of_pop + num_of_pop_per_split) / 2, sizeof(double));
    }

    pop_rcvbuf_split_weight = (double *)malloc(sizeof(double) * num_sendparams);
    pop_rcvbuf_split_delay = (double *)malloc(sizeof(double) * num_sendparams);

    pop_split_whole = (double *)malloc(sizeof(double) * num_of_pop_per_split * dimension);
    
    pop_sendbuf_nrn_weight = (double *)calloc(offset + num_sendparams, sizeof(double));
    pop_sendbuf_nrn_delay = (double *)calloc(offset + num_sendparams, sizeof(double));
    pop_rcvbuf_nrn_weight = (double *)malloc(num_of_weight_delay_per_procs * sizeof(double));
    pop_rcvbuf_nrn_delay= (double *)malloc(num_of_weight_delay_per_procs * sizeof(double));

    arFunvals_split_buf1 = (double *)calloc(num_of_pop_per_split, sizeof(double));
    arFunvals_split_buf2 = (double *)calloc(num_of_pop_per_split * num_of_procs_nrn + num_of_pop_per_split , sizeof(double));

    arFunvals_whole_buf = (double *)calloc(num_of_pop_per_split, sizeof(double));
    arFunvals_whole = (double *)calloc(num_of_pop + num_of_pop_per_split, sizeof(double));
    //x_temp = (double *)malloc(dimension * sizeof(double));
    x_temp = (double *)calloc(dimension, sizeof(double));

    sp_weight = (double *)malloc(sizeof(double) * mu);

    arFunval_rank = (int *)malloc(sizeof(int) * num_of_pop);
    sum_eachprocs = (double *)malloc(sizeof(double) * dimension);
    if(I_AM_ROOT_IN_MAIN){
	divider1 = 1.0;
    }else{
	divider1 = 1.0 / (num_of_pop_per_split * main_myid);
    }
    divider2 = 1.0 / (num_of_pop_per_split * (main_myid + 1));

    /*flgdiag off ver.*/
    dim_cov = (1 + dimension) * dimension / 2;
    /*flgdiag on ver.*/
    //dim_cov = dimension;
    sum_for_cov = (double *)malloc(sizeof(double) * dim_cov); 

    sum_reduce = (double *)malloc(sizeof(double) * dimension);
    sum_for_cov_reduce = (double *)malloc(sizeof(double) * dim_cov);
    
    initialX = (double *)malloc(sizeof(double) * dimension);
    if(initialX==NULL){ printf("memory allocation error for initialX \n"); return -1;}
    initialSigma = (double *)malloc(sizeof(double) * dimension);
    if(initialSigma==NULL){ printf("memory allocation error for initialSigma \n"); return -1;}
    restartSigma = (double *)malloc( sizeof(double) * dimension);
    //if(restartSigma == NULL){ printf("memory allocation error for restartSigma\n");}
    srand((unsigned)time(NULL));
    util = 1.0 / ((double)RAND_MAX + 2.0);
    for(i=0; i<dimension; ++i){
	initialX[i] = 3.0 + (7.0 - 3.0) * (double)(rand() + 1.0) * util;
	initialSigma[i] = (10.0 - 0.0) * 0.25;
	restartSigma[i] = restartSigma_defaults;
    }

    if(I_AM_ROOT_IN_MAIN){
	arFunvals = cmaes_init(&evo, dimension, initialX, initialSigma, seed, num_of_pop, mu, max_eval, max_iter, initfile);
	rand_box = evo.rand;
	for(i=0;i<mu;++i){
	    sp_weight[i] = evo.sp.weights[i];
	}
    }else{
	arFunvals = (double *)malloc(num_of_pop * sizeof(double));
	rand_box.rgrand = (long *)calloc(32, sizeof(long));
    }
    
    max_iter = (int)cmaes_Get(&evo, "MaxIter");
    max_eval = (int)cmaes_Get(&evo, "maxeval");

    /* setting MPI_Comm_split section*/
    key = 0;
    color = main_myid;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &splitcomm);
    MPI_Comm_rank(splitcomm, &split_myid);
    MPI_Comm_size(splitcomm, &split_size);

    MPI_Comm_dup(MPI_COMM_WORLD, &firstTimeWorld);


    MPI_Bcast(&rand_box.startseed, 1, MPI_LONG, root_process_main, firstTimeWorld);
    MPI_Bcast(&rand_box.aktseed, 1, MPI_LONG, root_process_main, firstTimeWorld);
    MPI_Bcast(&rand_box.aktrand, 1, MPI_LONG, root_process_main, firstTimeWorld);
    MPI_Bcast(rand_box.rgrand, 32, MPI_LONG, root_process_main, firstTimeWorld);
    MPI_Bcast(&rand_box.flgstored, 1, MPI_SHORT, root_process_main, firstTimeWorld);
    MPI_Bcast(&rand_box.hold, 1, MPI_DOUBLE, root_process_main, firstTimeWorld);
    rand_box.startseed += main_myid;
    rand_box.aktseed += main_myid;
    MPI_Bcast(sp_weight, mu, MPI_DOUBLE, root_process_main, firstTimeWorld);
    
    /*setting the argv for spawn*/
    sprintf(neuron_argv[3],"COLOR=%d",color);
    if(neuron_argv_size > 6){
    	sprintf(neuron_argv[5], "NCELL_CMAES=%d", network_size);
    	neuron_argv[neuron_argv_size-1] = NULL;
    }else{
    	neuron_argv[5] = NULL;
    }
    /*********caution********/
    /*probably the number of split comm is equal to the number of the main process, so the if sequence below may be unnecessary.*/
    /*********caution********/
    /*execute MPI_Comm_spawn (from below, until sending finalize signal to nrncomm, the program should not end.*/
    if( I_AM_ROOT_IN_SPLIT){
	MPI_Comm_spawn(specials, neuron_argv, num_of_procs_nrn, MPI_INFO_NULL, 0, splitcomm, &intercomm, MPI_ERRCODES_IGNORE);
	MPI_Intercomm_merge(intercomm, 0, &nrn_comm);
	MPI_Comm_size(nrn_comm, &spawn_size);
	MPI_Comm_rank(nrn_comm, &spawn_myid);
    }
    /*if you want to send information to NEURON in setting, send here*/
    if( I_AM_ROOT_IN_SPLIT){
	send_count = 1;//3;/*template*/
	double info[] = {1,0, 1.0, 1.0};/*template*/
	info[0] = num_of_pop_per_split;
	//info[1] and info[2] are setted default now..( may change in future plan)
       	MPI_Bcast_to_NEURON(info, send_count, MPI_DOUBLE, root_process_spawn, nrn_comm);
    }
    fflush(stdout);

    t_start = MPI_Wtime();
    printf("main_myid = %d \n", main_myid);

    /*Start main section of estimation*/
    while(1){
	if( I_AM_ROOT_IN_MAIN){
	    rgD = (double *)cmaes_GetPtr(&evo, "diag(C)");
	    x_mean = (double *)cmaes_GetPtr(&evo, "xmean");
	    cmaes_SamplePopulation_diag_dist_update(&evo);
	    sigma = (double)cmaes_Get(&evo, "sigma");
	}else{
	    rgD = (double *)malloc(dimension * sizeof(double));
	    x_mean = (double *)malloc(dimension * sizeof(double));
	}

	MPI_Bcast(rgD, dimension, MPI_DOUBLE, root_process_main, firstTimeWorld);
	MPI_Bcast(x_mean, dimension, MPI_DOUBLE, root_process_main, firstTimeWorld);
	MPI_Bcast(&sigma, 1, MPI_DOUBLE, root_process_main, firstTimeWorld);

	sigmasquare_div = 1.0 / (sigma * sigma);
	
	for(i = 0; i < num_of_pop_per_split; ++i){
	    pop_dist = (double *)cmaes_SamplePopulation_diag_dist(rgD, sigma, x_mean, dimension, rand_box);
	    for(j=0;j<dimension;++j){
		pop_split_whole[i*dimension + j] = pop_dist[j];
	    }
	    my_boundary_transformation(&my_boundaries, pop_dist, x_temp, main_myid);
	    for(j = 0; j< (dimension / 2); ++j){
		/* pop_split_whole[i * dimension + j] = x_temp[j]; */
		/* pop_split_whole[i * dimension + j + dimension / 2] = x_temp[j + dimension / 2]; */
		pop_rcvbuf_split_weight[i * dimension / 2 + j] = x_temp[j];
		pop_rcvbuf_split_delay[i * dimension / 2 + j] = x_temp[j + dimension / 2];
	    }
	}
	/*past not distributed ver.*/
	
    	/* if( I_AM_ROOT_IN_MAIN ){ */
    	/*     pop = cmaes_SamplePopulation(&evo);/\*do not change content of pop*\/ */
	/*     printf("x_temp\'s address is %p (myid = %d)\n", x_temp, main_myid); */
    	/*     for(i=0;i<num_of_pop;++i){ */
    	/* 	my_boundary_transformation(&my_boundaries, pop[i], x_temp, main_myid); */
    	/* 	for(j=0;j<dimension;++j){ */
    	/* 	    pop_sendbuf_split_whole[i * dimension + j] = x_temp[j]; */
    	/* 	} */
    	/*     } */
	/*     for(i=0; i<num_of_pop; ++i){ */
	/* 	for(j=0; j<(dimension/2); ++j){ */
	/* 	    pop_sendbuf_split_weight[i * dimension / 2 + j] = pop_sendbuf_split_whole[i * dimension + j]; */
	/* 	    pop_sendbuf_split_delay[i * dimension / 2 + j] = pop_sendbuf_split_whole[i * dimension + j + dimension / 2]; */
	/* 	} */
	/*     } */
	/*     printf("end of generate pops\n");	 */
	/* } */
	/* //when you spawn, it can be that scatter in MPI_COMM_WORLD (but the line below is temporaly) */
	/* MPI_Scatter(pop_sendbuf_split_weight, num_sendparams, MPI_DOUBLE, pop_rcvbuf_split_weight, num_sendparams, MPI_DOUBLE, root_process_main, firstTimeWorld); */
	/* MPI_Scatter(pop_sendbuf_split_delay, num_sendparams, MPI_DOUBLE, pop_rcvbuf_split_delay, num_sendparams, MPI_DOUBLE, root_process_main, firstTimeWorld); */
	/*MPI_Barrier(MPI_COMM_WORLD);*/
	if(I_AM_ROOT_IN_SPLIT){
	    for(k=0;k<num_of_procs_nrn;k++){
		for(i=0;i<num_of_pop_per_split;i++){
		    for(j=0;j<num_of_one_gene_weight_or_delay;j++){
			pop_sendbuf_nrn_weight[offset + (num_of_pop_per_split * num_of_one_gene_weight_or_delay) * k + (num_of_one_gene_weight_or_delay) * i + j]= pop_rcvbuf_split_weight[ (dimension / 2) * i + j + num_of_one_gene_weight_or_delay * k];
			pop_sendbuf_nrn_delay[offset + (num_of_pop_per_split * num_of_one_gene_weight_or_delay) * k + (num_of_one_gene_weight_or_delay) * i + j]= pop_rcvbuf_split_delay[ (dimension / 2) * i + j + num_of_one_gene_weight_or_delay * k];
		    }
		}
	    }
    	    /*evaluate the new searching points*/
	    /*fatal error in this scatter section*/
    	    MPI_Scatter(pop_sendbuf_nrn_weight, num_of_weight_delay_per_procs, MPI_DOUBLE, pop_rcvbuf_nrn_weight, num_of_weight_delay_per_procs, MPI_DOUBLE, root_process_spawn, nrn_comm);
	    MPI_Scatter(pop_sendbuf_nrn_delay, num_of_weight_delay_per_procs, MPI_DOUBLE, pop_rcvbuf_nrn_delay, num_of_weight_delay_per_procs, MPI_DOUBLE, root_process_spawn, nrn_comm);

	    MPI_Gather(pop_split_whole, num_of_pop_per_split * dimension, MPI_DOUBLE, pop_sendbuf_split_whole, num_of_pop_per_split * dimension, MPI_DOUBLE, root_process_main, firstTimeWorld);
	    if(I_AM_ROOT_IN_MAIN){
		for(i = 0; i < num_of_pop; ++i){
		    for(j = 0; j < dimension; ++j){
			evo.rgrgx[i][j] = pop_sendbuf_split_whole[i * dimension + j];
			//printf("evo.rgrgx[%d][%d] = %lf\n", i, j, evo.rgrgx[i][j]);
		    }
		}
	    }
	    
	    
    	    /* wait for NEURON simulation in worker nodes */
    	    MPI_Gather(arFunvals_split_buf1, num_of_pop_per_split, MPI_DOUBLE, arFunvals_split_buf2, num_of_pop_per_split, MPI_DOUBLE, root_process_spawn, nrn_comm);
	    MPI_Barrier(firstTimeWorld);
	    printf("end of gather\n");
    	    for(i=0;i<num_of_pop_per_split; ++i){
    		arFunvals_whole_buf[i] = arFunvals_split_buf2[i + num_of_pop_per_split];
    	    }
	}
	/*then collect information of all split process and update (temporaly setting parameter)*/
	//MPI_Gather(arFunvals_whole_buf, num_of_pop_per_split, MPI_DOUBLE, arFunvals_whole, num_of_pop_per_split, MPI_DOUBLE, root_process_main, MPI_COMM_WORLD);
	
	MPI_Allgather(arFunvals_whole_buf, num_of_pop_per_split, MPI_DOUBLE, arFunvals_whole, num_of_pop_per_split, MPI_DOUBLE, MPI_COMM_WORLD);
	Sort_score_index(arFunvals_whole, arFunval_rank, num_of_pop);
	
	for(i=0; i<dimension; ++i){
	    sum_eachprocs[i] = 0;
	}

	for(i=0; i<dim_cov; ++i){
	    sum_for_cov[i] = 0;
	}

	/*calculate for updateDistribution*/
	for(i = 0; i < mu; ++i){
	    for(j = 0; j < dimension; ++j){
		sum_eachprocs[j] += ((arFunval_rank[i]!=0)*((int)(arFunval_rank[i] * divider1)) * ((int)(arFunval_rank[i] * divider2)==0) + (arFunval_rank[i]==0) * (I_AM_ROOT_IN_MAIN)) * sp_weight[i] * pop_split_whole[arFunval_rank[i]%num_of_pop_per_split * dimension + j];
	    }
	}

	/*calculate for Adapt_C2, flg_diag_off_ver.*/
	for(k = 0; k < mu; ++k){
	    sum_cov_counter=0;
	    for(i = 0; i < dimension; ++i){
		for(j = 0; j <= i; ++j){
		    sum_for_cov[sum_cov_counter] += ((arFunval_rank[k]!=0) * ((int)(arFunval_rank[k] * divider1)) * ((int)(arFunval_rank[k] * divider2)==0) + (arFunval_rank[k]==0) * (I_AM_ROOT_IN_MAIN) )* sp_weight[k] * (pop_split_whole[arFunval_rank[k]%num_of_pop_per_split * dimension + i] - x_mean[i]) *  (pop_split_whole[arFunval_rank[k]%num_of_pop_per_split * dimension + j] - x_mean[j]) * sigmasquare_div;
		    ++sum_cov_counter;
								     
		}
	    }
	}

	/*send calc results to root_process_main*/
	MPI_Reduce(sum_eachprocs, sum_reduce, dimension, MPI_DOUBLE, MPI_SUM, root_process_main, firstTimeWorld);
	MPI_Reduce(sum_for_cov, sum_for_cov_reduce, dim_cov, MPI_DOUBLE, MPI_SUM, root_process_main, firstTimeWorld);
	
	
	if(I_AM_ROOT_IN_MAIN){
	    for(i=0;i<num_of_pop;++i){
		arFunvals[i] = arFunvals_whole[i];
		if(arFunvals[i] == 0){
		    printf("communication error\n");
		    break;
		}
	    }
    	    /*update the search distribution used for cmaes_sampleDistribution()*/
    	    //cmaes_UpdateDistribution(&evo, arFunvals); /*assume that pop[] has not been modified*/
	    cmaes_UpdateDistribution_dist(&evo, arFunvals, sum_reduce, sum_for_cov_reduce);
    	}

	/* if(I_AM_ROOT_IN_MAIN){ */
	/*     for(i=0;i<num_of_pop;++i){ */
	/* 	for(j=0;j<dimension;++j){ */
	/* 	    printf("(past)evo.rgrgx[%d][%d] = %lf\n", i, j, evo.rgrgx[i][j]); */
	/* 	} */
	/*     } */
	/* } */
	
    	fflush(stdout);
    	/*termination*/
    	if(I_AM_ROOT_IN_MAIN){
    	    if(cmaes_TestForTermination(&evo)){
    		flg_termination = 1;
    	    }
    	    send_count = 1;
	}
	/* broadcast flg to MAIN PROCS*/
	MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, root_process_main, firstTimeWorld);
		    
	if(I_AM_ROOT_IN_SPLIT){
    	    MPI_Bcast_to_NEURON(&flg_termination, 1, MPI_DOUBLE, root_process_split, nrn_comm);
    	    //MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, root_process_split, splitcomm);
    	    if((int)flg_termination){
    		break;
    	    }
    	}/*cmaes termination*/

	/*for restart */
	if((int)cmaes_Get(&evo, "generation") == gen_restart[n_run] && I_AM_ROOT_IN_MAIN){
	    generation = cmaes_Get(&evo, "generation");
	    countevals = cmaes_Get(&evo, "eval");
	    restartX = (double *)cmaes_GetPtr(&evo, "xmean");
	    /* for(k=0;k<dimension;k++){ */
	    /* 	restartSigma[k] = cmaes_Get(&evo, "sigma"); */
	    /* 	printf("restartSigma[%d] = %lf\n",k, restartSigma[k]); */
	    /* }*/
	    //x_temp = (double *)cmaes_GetPtr(&evo, "xbestever");
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
    
    MPI_Barrier(MPI_COMM_WORLD);
    t_end = MPI_Wtime();
    if( I_AM_ROOT_IN_MAIN){
    	printf("#Stop:\n%s\n", cmaes_TestForTermination(&evo)); /*print termination reason*/
    	printf("\n# operation finished.\n# elapsed time: %f\n #fbest: %f\n #xbest:", t_end - t_start, cmaes_Get(&evo, "fbestever"));
    	my_boundary_transformation(&my_boundaries, (double *)cmaes_GetPtr(&evo, "xbestever"), x_temp, main_myid);
    	printGene(stdout, x_temp, dimension);
	printf("\n");
	fflush(stdout);
    }

    /* /\* for communication test between nrncomms*\/ */
    /* test_sendbuf = (double *)calloc(20, sizeof(double)); */
    /* test_rcvbuf = (double *)malloc(sizeof(double) * 4); */
    /* test_arFunval1 = (double *)calloc(4, sizeof(double)); */
    /* test_arFunval2 = (double *)calloc(20, sizeof(double)); */
    /* while(1){ */
    /* 	if(I_AM_ROOT_IN_SPLIT){ */
    /* 	    for(i=0;i<20;i++){ */
    /* 		test_sendbuf[i] = i + 10; */
    /* 	    } */
    /* 	    MPI_Scatter(test_sendbuf, 4, MPI_DOUBLE, test_rcvbuf, 4, MPI_DOUBLE, root_process_spawn, nrn_comm); */

    /* 	    MPI_Gather(test_arFunval1, 4, MPI_DOUBLE, test_arFunval2, 4, MPI_DOUBLE, root_process_spawn, nrn_comm); */
    /* 	    for(i=0;i<20;i++){ */
    /* 		printf("test_arFunval2[%d] = %lf\n",i, test_arFunval2[i]); */
    /* 	    } */
    /* 	} */
    /* 	if(I_AM_ROOT_IN_SPLIT){ */
    /* 	    flg_termination = 1; */
    /* 	    send_count = 1; */
    /* 	    MPI_Bcast_to_NEURON(&flg_termination, 1, MPI_DOUBLE, root_process_split, nrn_comm);/\*before splitcomm -> after nrn_comm*\/ */
    /* 	} */
    /* 	printf("start MPI_Bcast\n"); */
    /* 	MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, root_process_split, splitcomm); */
    /* 	if((int)flg_termination){ */
    /* 	    break; */
    /* 	} */
    /* }/\*end of communication test loop*\/ */
    /* /\*free the allocate memory for test loop*\/ */
    /* free(test_sendbuf); */
    /* free(test_rcvbuf); */
    /* free(test_arFunval1); */
    /* free(test_arFunval2); */
    /* /\*end of free test memory*\/ */

    printf("end of loop\n");
    fflush(stdout);
    
    MPI_Barrier(MPI_COMM_WORLD);
    

    printf("exec loop time is %lf\n", t_end - t_start);
    
    /* finalize the process (free the memory)*/
    if(I_AM_ROOT_IN_MAIN){
	cmaes_exit(&evo);
    }
    my_boundary_transformation_exit(&my_boundaries);//after define loadRangefile, uncomment!!
    free(x_temp);
    free(initialX);
    free(initialSigma);
    if(I_AM_ROOT_IN_MAIN){
	free(pop_sendbuf_split_whole);
	//free(pop_sendbuf_split_delay);
	//free(pop_sendbuf_split_weight);
    }
    free(pop_sendbuf_nrn_delay);
    free(pop_sendbuf_nrn_weight);
    free(pop_rcvbuf_split_delay);
    free(pop_rcvbuf_split_weight);
    free(pop_rcvbuf_nrn_weight);
    free(pop_rcvbuf_nrn_delay);
    free(arFunvals_split_buf1);
    free(arFunvals_split_buf2);
    free(arFunvals_whole_buf);
    free(arFunvals_whole);

    for(i=0;i<6;i++){
	//free(neuron_argv[i]);
    }
    //free(neuron_argv);

    printf("end of the free memory section (myid = %d)\n", main_myid);
    
    MPI_Finalize();

    return 0;
}
