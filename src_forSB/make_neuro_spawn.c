#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>


/*definition of the wrapper of MPI_Bcast (function for communication with NEURON) */
void MPI_Bcast_to_NEURON(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    /* datatype: now support MPI_DOUBLE only*/
    MPI_Bcast(&count, 1, MPI_INT, root, comm);
    MPI_Bcast(buffer, count, datatype, root, comm);
    fflush(stdout);
}/*MPI_Bcast_to_NEURON*/


int main(int argc, char **argv){
    int i, j, k;
    MPI_Comm parentcomm, spawn_parent_comm, intercomm, nrn_comm;
    int main_size, main_myid;
    int spawn_parent_size, spawn_parent_rank;
    int spawn_size, spawn_myid;
    int num_of_procs_nrn = 2; // change in your environment (local, K, Fx100)
    int dimension_per_one_nrnproc, num_of_only_weights,  num_of_params_of_one_nrnrproc;
    int dimension_per_one_nrnproc_w;
    int num_of_weights_per_one_nrnproc;
    int num_of_my_pop, dimension;
    int num_sendparams_from_parent, num_sendparams_to_NEURON;
    int offset;

    double *pop_sendbuf_whole, *pop_rcvbuf_whole;
    double *pop_sendbuf_child_weight;
    double *pop_rcvbuf_child_weight;

    double *pop_sendbuf_nrn_weight;
    double *pop_rcvbuf_nrn_weight;

    double *arFunvals_child_buf1, *arFunvals_child_buf2;
    double *arFunvals_whole_buf, *arFunvals_whole;
    
    int send_count;
    double flg_termination;
    
    char specials[] = "../hocfile_forSB/x86_64/special";
    char *neuron_argv[] = {"-mpi", "-nobanner", "../hocfile_forSB/networkSimulation.hoc", NULL};

    char connection_data[] = "../data/conMat.txt";
    FILE *fp;
    int dim_conMat=2; //ncell in NEURON
    int **conMat;
    int num_of_cell_combination = dim_conMat * dim_conMat;
    int gene_id;

    /* initialize the MPI process*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    /* make new intracommunicator for communicating with the parent process*/
    MPI_Comm_get_parent(&parentcomm);
    MPI_Intercomm_merge(parentcomm, 1, &spawn_parent_comm);
    MPI_Comm_size(spawn_parent_comm, &spawn_parent_size);
    MPI_Comm_rank(spawn_parent_comm, &spawn_parent_rank);

    /* recieve the information of the scale of population from parent process*/
    /* required: num_of_my_pop, dimension, */
    if(argc>2){
	num_of_my_pop = atoi(argv[1]);
	dimension = atoi(argv[2]);
    }else{
	num_of_my_pop = 32;
	dimension = 72;
    }

    /* set variables for communication */
    num_sendparams_from_parent = dimension * num_of_my_pop;/* separate weight and delay ver. (for non-separate, delete / 2*/

    dimension_per_one_nrnproc = num_of_cell_combination / num_of_procs_nrn;
    dimension_per_one_nrnproc_w = dimension_per_one_nrnproc;
    num_of_only_weights = num_of_cell_combination * num_of_my_pop;
    num_of_weights_per_one_nrnproc = num_of_only_weights / num_of_procs_nrn;
    offset = num_of_weights_per_one_nrnproc;
    num_sendparams_to_NEURON = num_of_cell_combination * num_of_my_pop;

    /* allocate communication arrays (some arrays are not allocated because they do not need in these processes) */
    

    pop_rcvbuf_whole = (double *)malloc(sizeof(double) * num_sendparams_from_parent);

    pop_sendbuf_nrn_weight = (double *)malloc(sizeof(double) * (offset + num_sendparams_to_NEURON));
    pop_rcvbuf_nrn_weight = (double *)malloc(sizeof(double) * num_of_weights_per_one_nrnproc);

    arFunvals_child_buf1 = (double *)calloc(num_of_my_pop, sizeof(double));
    arFunvals_child_buf2 = (double *)calloc(num_of_my_pop * num_of_procs_nrn + num_of_my_pop, sizeof(double));

    arFunvals_whole_buf = (double *)calloc(num_of_my_pop, sizeof(double));
    arFunvals_whole = (double *)calloc(num_of_my_pop * spawn_parent_size, sizeof(double));

    conMat = (int **)malloc(sizeof(int *) * dim_conMat);
    for(i=0;i<dim_conMat;++i){
	conMat[i] = (int *)malloc(sizeof(int) * dim_conMat);
    }
    if((fp=fopen(connection_data, "r"))==NULL){
	printf("connection file open error\n");
	exit(EXIT_FAILURE);
    }
    for(i=0;i<dim_conMat;++i){
	for(j=0;j<dim_conMat;++j){
	    fscanf(fp, "%d", &conMat[i][j]);
	}
    }
    fclose(fp);

    /* make spawn of NEURON procs and make new intracommunicator 'nrn_comm'*/
    /* when it does not work, uncomment the below sentense*/
    /* for(i=0; i < 8; ++i){ */
    /* 	if(spawn_parent_rank%8 == i){ */
    MPI_Comm_spawn(specials, neuron_argv, num_of_procs_nrn, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(intercomm, 0, &nrn_comm);
    MPI_Comm_size(nrn_comm, &spawn_size);
    MPI_Comm_rank(nrn_comm, &spawn_myid);
    /* 	} */
    /* 	MPI_Barrier(MPI_COMM_WORLD); */
    /* } */

    /* send information of the scale of population to NEURON process*/
    send_count = 1;
    double info[] = {1.0, 1.0, 1.0};
    info[0] = num_of_my_pop; 
    MPI_Bcast_to_NEURON(info, send_count, MPI_DOUBLE, 0, nrn_comm);

    /* infinite loop for estimation (receive information of genes and pass it to NEURON process*/
    while(1){
	/* recieve the gene information from parent process*/
	MPI_Scatter(pop_sendbuf_whole, num_sendparams_from_parent, MPI_DOUBLE, pop_rcvbuf_whole, num_sendparams_from_parent, MPI_DOUBLE, 0, spawn_parent_comm);
	/* transfrom the data structure to suitable manner for communication */
	
	gene_id = 0;
	for(k=0;k<num_of_my_pop;++k){
	    for(i=0;i<dim_conMat;++i){
		for(j=0;j<dim_conMat;++j){
		    pop_sendbuf_nrn_weight[offset + num_of_cell_combination * k + dimension_per_one_nrnproc * i + j] = pop_rcvbuf_whole[gene_id] * conMat[i][j];
		    gene_id += conMat[i][j];
		}
	    }
	}

	/* for(k=0;k<num_of_my_pop;++k){ */
	/*     for(i=0;i<dim_conMat;++i){ */
	/* 	for(j=0;j<dim_conMat;++j){ */
	/* 	    printf("w[%d][%d] = %lf \t", i, j, pop_sendbuf_nrn_weight[offset + num_of_cell_combination * k + dimension_per_one_nrnproc * i + j]); */

	/* 	} */
	/* 	printf("\n"); */
	/*     } */
	/*     printf("\n\n"); */
	/* } */
	

	MPI_Barrier(MPI_COMM_WORLD);

	/*pass the gene information to NEURON process*/
	MPI_Scatter(pop_sendbuf_nrn_weight, num_of_weights_per_one_nrnproc, MPI_DOUBLE, pop_rcvbuf_nrn_weight, num_of_weights_per_one_nrnproc, MPI_DOUBLE, 0, nrn_comm);

	/* wait for finish calculation in NEURON*/
	MPI_Gather(arFunvals_child_buf1, num_of_my_pop, MPI_DOUBLE, arFunvals_child_buf2, num_of_my_pop, MPI_DOUBLE, 0, nrn_comm);

	/* tally the score information of my process */
	for(i=0; i<num_of_my_pop; ++i){
	    arFunvals_whole_buf[i] = arFunvals_child_buf2[i + num_of_my_pop];
	    //printf("arFunvals_whole_buf[%d] = %lf\n", i, arFunvals_whole_buf[i]);
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	/* pass the score information to parent process*/
	MPI_Gather(arFunvals_whole_buf, num_of_my_pop, MPI_DOUBLE, arFunvals_whole, num_of_my_pop, MPI_DOUBLE, 0, spawn_parent_comm);
	/* recieve and pass the information whether the loop is terminated or not */
	MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, 0, spawn_parent_comm);
	MPI_Bcast_to_NEURON(&flg_termination, 1, MPI_DOUBLE, 0, nrn_comm);
	/* termination*/
	if((int)flg_termination){
	    break;
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* free the allocated memory */
    free(pop_rcvbuf_child_weight);
    free(pop_sendbuf_nrn_weight);
    free(pop_rcvbuf_nrn_weight);
    free(arFunvals_child_buf1);
    free(arFunvals_child_buf2);
    free(arFunvals_whole_buf);
    free(arFunvals_whole);
   
    for(i=0;i<dim_conMat;++i){
	free(conMat[i]);
    }
    free(conMat);
 
    MPI_Finalize();

    return 0;
}
