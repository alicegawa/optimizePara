#include "mplmcma.h"

int main(int argc, char **argv){
    int i;
    int main_myid, main_size;
    char command[] = "./make_neuro_spawn";
    char **spawn_argv;
    int spawn_argv_size=7;
    int num_of_pop_per_child;
    int spawn_numprocs = 7;/*the number of NEURON circuits*/
    MPI_Comm spawn_comm, intercomm, parentcomm;
    int spawn_size, spawn_myid;
    int num_of_nrn_procs=4;
    char *exec_prog=NULL;
    int dim_con_mat;
    char *connection_data=NULL;

    int N;                  //      problem dimension
    int lambda;
    int mu;
    double ccov;// = 1/(10*log((double)N+1.0));//      learning rate for covariance matrix, e.g., 1/(10*log(N+1))
    double xmin = -5;//     x parameters lower bound
    double xmax = 5;//      x parameters upper bound
    double *xmin_vec, *xmax_vec;
    double *xmin_tmp, *xmax_tmp;
    int nvectors;  //      number of stored direction vectors, e.g., nvectors = 4+floor(3*log(N))
    int maxsteps;        //      target number of generations between vectors, e.g., maxsteps = nvectors
    double cc;     // learning rate for mean vector's evolution path, e.g., cc = 1/nvectors
    double val_target = 0.25;       // target success rate for new population, e.g., 0.25
    double sigma = 0.5*(xmax - xmin);       // initial step-size, e.g., 0.5
    double *sigma_vec;
    double c_s = 0.3;       //      decay factor for step-size adaptation, e.g., 0.3
    double target_f = 1e-4;        // target fitness function value, e.g., 1e-10
    int maxevals = 5;            // maximum number of function evaluations allowed, e.g., 1e+6
    int inseed = 1;         // initial seed for random number generator, e.g., 1
    int algorithmType = 10; // 0:LMCMA, 10:LMCMA_fixed+features, 1:sepCMA, 2:Cholesky-CMA, 3: baseline CMA-ES
    int printToFile = 0; // 1 or 0

    int sample_symmetry = 0; // use 1 to make the algorithm 2 times faster (in CPU time) and perform better (in terms of number of evaluations)

    double output[2];

    char range_filename[256] = "../data/params_onlyWeight.txt";
    mplmcma_t t;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &main_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

    MPI_Comm_get_parent(&parentcomm);

    xmin_tmp = (double *)malloc(sizeof(double) * MAX_NUM_PARAM);
    xmax_tmp = (double *)malloc(sizeof(double) * MAX_NUM_PARAM);


    
    spawn_argv = (char **)malloc(sizeof(char *) * spawn_argv_size);
    for(i=0;i<spawn_argv_size;++i){
      spawn_argv[i] = (char *)malloc(sizeof(char) * 256);
    }
    spawn_argv[spawn_argv_size-1] = NULL;
    exec_prog = (char *)malloc(sizeof(char) * 256);
    if(exec_prog==NULL){
      printf("memory allocation error occurs @{exec_prog}\n");
    }
    connection_data = (char *)malloc(sizeof(char) * 256);
    if(connection_data==NULL){
      printf("memory allocation error occurs @{connection_data}\n");
    }

    /* for test input*/
    if(argc > 7){
      dim_con_mat = atoi(argv[7]);
      sprintf(connection_data, "%s", argv[8]);
      sprintf(range_filename, "%s", argv[9]);
    }
    
    loadRangeFile(range_filename, xmin_tmp, xmax_tmp, &t.N);
    N = t.N;

    if((dim_con_mat * dim_con_mat)!=N){
      dim_con_mat = 1;
      while(1){
	if((dim_con_mat * dim_con_mat)>=N){
	  break;
	}else{
	  dim_con_mat++;
	}
      }
    }
    printf("dim_con_mat = %d\n", dim_con_mat);
    
    if(argc < 2){
      lambda = 4+floor(3*log((double)N)); //      population size, e.g., 4+floor(3*log(N));
      mu = (int)(lambda / 2); // number of parents, e.g., floor(lambda /2);
    }else{
      t.lambda=lambda = atoi(argv[1]);
      t.mu = mu = atoi(argv[2]);
      t.spawn_numprocs = spawn_numprocs = atoi(argv[3]);
      t.maxevals = maxevals = atoi(argv[4]);
      maxevals *= lambda;
      t.maxevals *= lambda;
      t.num_of_nrn_procs = atoi(argv[5]);
      sprintf(exec_prog, "%s", argv[6]);
    }
    
    t.num_of_pop_per_child = num_of_pop_per_child = lambda / spawn_numprocs;
    t.nvectors = nvectors = lambda;
    t.maxsteps = maxsteps = nvectors;
    t.cc = cc = (double)(1.0/(double)nvectors);
    t.ccov = ccov = 1/(10*log((double)N+1.0));
    maxevals *= lambda;
    t.maxevals *= lambda;
    printf("N = %d, lambda = %d, mu = %d\n", N, lambda, mu);
    sprintf(spawn_argv[0], "%d", num_of_pop_per_child);
    sprintf(spawn_argv[1], "%d", N);
    sprintf(spawn_argv[2], "%d", num_of_nrn_procs);
    sprintf(spawn_argv[3], "%s", exec_prog);
    sprintf(spawn_argv[4], "%d", dim_con_mat);
    sprintf(spawn_argv[5], "%s", connection_data);

    xmin_vec = (double *)malloc(sizeof(double) * N);
    xmax_vec = (double *)malloc(sizeof(double) * N);

    init_hyperparam_arrays(&t);

    for(i=0;i<N;++i){
      t.xmin[i] = xmin_vec[i] = xmin_tmp[i];
      t.xmax[i] = xmax_vec[i] = xmax_tmp[i];
    }
    free(xmin_tmp); xmin_tmp = NULL;
    free(xmax_tmp); xmax_tmp = NULL;

    sigma_vec = (double *)malloc(sizeof(double) * N);
    for(i=0;i<N;++i){
	t.sigma[i] = sigma_vec[i] = 0.5 * (xmax_vec[i] - xmin_vec[i]);
    }

    printf("start comm spawn\n");
    MPI_Comm_spawn(command, spawn_argv, spawn_numprocs, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(intercomm, 0, &spawn_comm);
    MPI_Comm_size(spawn_comm, &spawn_size);
    MPI_Comm_rank(spawn_comm, &spawn_myid);
    
    t.spawn_comm = spawn_comm;
    t.printToFile = 0;
    t.sample_symmetry = 0;
    t.inseed = inseed;
    t.val_target = val_target;
    t.c_s = c_s;
    t.target_f = 1e-5;

    printf("start LMCMA\n");
    LMCMA(&t);
    //LMCMA(t.N, t.lambda, t.mu, t.ccov, t.xmin, t.xmax, t.nvectors, t.maxsteps, t.cc, t.val_target, t.sigma, t.c_s, t.target_f, t.maxevals, t.inseed, output, t.printToFile, t.sample_symmetry, t.spawn_comm, t.spawn_numprocs, &t);
    //LMCMA(N, lambda, mu, ccov, xmin_vec, xmax_vec, nvectors, maxsteps, cc, val_target, sigma_vec, c_s, target_f, maxevals, inseed, output, printToFile, sample_symmetry, spawn_comm, spawn_numprocs);
    printf("end of LMCMA\n");

    MPI_Comm_free(&intercomm);
    MPI_Comm_free(&spawn_comm);

    free(xmax_vec); xmax_vec = NULL;
    free(xmin_vec); xmin_vec = NULL;
    free(exec_prog); exec_prog = NULL;
    free(connection_data); connection_data=NULL;
    free_hyperparam_arrays(&t);

    MPI_Finalize();
    return 0;
}
