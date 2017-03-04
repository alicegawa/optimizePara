#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "my_boundary_transformation.h"

#define TEXT_BUFFER_SIZE 2048
#define MAX_NUM_PARAM 10000
#define MAX_NUM_TARGET 5

#define USE_ZIGGURAT 0

typedef struct 
/* random_t 
 * sets up a pseudo random number generator instance 
 */
{
  /* Variables for Uniform() */
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int *rgrand;
  
  /* Variables for Gauss() */
  short flgstored;
  double hold;
} random_t; 

typedef struct 
{
  double value;
  int id;
} sortedvals;

typedef struct
{
  random_t ttime;
  double*	func_tempdata;
  double*	x_tempdata;
  double*	rotmatrix;
  double* func_shiftxi;
  time_t	time_tic_t;
  clock_t	time_tic_c;
  time_t	time_toc_t;
  clock_t	time_toc_c;
} global_t;

typedef struct
{
  /*hyper params*/
  int N;
  int lambda;
  int mu;
  double ccov;
  double *xmin, *xmax;
  int nvectors;
  int maxsteps;
  double cc;
  double val_target;
  double *sigma;
  double c_s;
  double target_f;
  int maxevals;
  int inseed;
  double output[2];
  int printToFile;
  int sample_symmetry;
  MPI_Comm spawn_comm;
  int num_of_spawn_comm;
  int num_of_pop_per_spawn;
  int num_of_nrn_procs;
  double flg_termination;
  int num_of_pop_per_child;
  int spawn_numprocs;

  /*lmcma arrays*/
  double *arx;
  double *v_arr;
  double *pc_arr;
  double *pc;
  double *xmean;
  double *xold;
  double* xbestever;
  double* z;
  double* Az;
  double* Av;

  double* weights;
  int* iterator;
  double* arfitness;
  double* prev_arfitness;
  int* arindex;
  double* mixed;
  int* ranks;
  int* ranks_tmp;
  double* Nj_arr;
  double* Lj_arr;
  sortedvals* arr_tmp;
  int* ti;
  int* vec;
   
  /*for communicate with the child process*/
  double *scatter_sendvec_w, *scatter_sendvec_d;
  double *scatter_sendvec;
  double *scatter_rcvvec_w, *scatter_rcvvec_d;
  double *scatter_rcvvec;
  double *arFunvals_send, *arFunvals_rcv;

  double *x_temp;
  unsigned int *flg_log;

}mplmcma_t;

/* proto-type declaration*/
void init_gt(global_t* gt);
void free_gt(global_t* gt);
void random_exit(random_t *t);
long random_Start(random_t *t, long unsigned inseed);
long random_init(random_t *t, long unsigned inseed);
double random_Uniform(random_t *t);
double random_Gauss(random_t *t);
void time_tic(global_t *t);
double time_tictoc(global_t *t);
double time_toc(global_t *t);

void matrix_eye(double *M, int m);
void matrix_mult_vector(double *res, double *a, double *b, int m);
void matrix_mult_matrix(double *res, double *a, double *b, int m);
void vector_mult_vector(double *res, double *a, double *b, int m);
void vector_mult_matrix(double *res, double *a, double *b, int m);
double vector_prod(double *a, double *b, int m);
void generateRotationMatrix(double *B, int N, double *tmp1, random_t* rnd);

double minv(double a, double b);
double maxv(double a, double b);

double fsphere(double* x, int N);
double felli(double* x, int N);
double felli_fast(double* x, int N, global_t* t);
double fdiscus(double* x, int N);
double fcigar(double* x, int N);
void getRotatedX(double* x, int N, global_t* t);
double frosen(double* x, int N);

int compare(const void * a, const void * b);
void myqsort(int sz, double* arfitness, int* arindex, sortedvals* arr);

void invAz(int N, double* Av, int iterator_sz, int* iterator, double* v_arr, double* Lj_arr, double K);

int printGene(FILE *fp, const double *x, int dimension);
int printGene2(FILE *fp, const double *x, int dimension);

/* void LMCMA(int N, int lambda, int mu, double ccov, double *xmin, double *xmax, int nvectors, */
/* 	   int maxsteps, double cc, double val_target, double *sigma, double c_s, double target_f,  */
/* 	   int maxevals, int inseed, double* output, int printToFile, int sample_symmetry, MPI_Comm spawn_comm, int num_of_spawn_comm); */
int loadRangeFile(char *filename, double *xmin_vec, double *xmax_vec, int *N);

void nullset2arrays(mplmcma_t *t);
int init_lmcma_arrays(mplmcma_t *t);
void free_lmcma_arrays(mplmcma_t *t);

int init_hyperparam_arrays(mplmcma_t *t);
void free_hyperparam_arrays(mplmcma_t *t);
