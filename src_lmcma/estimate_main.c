#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "my_boundary_transformation.h"

#define TEXT_BUFFER_SIZE 2048
#define MAX_NUM_PARAM 150
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

void init_gt(global_t* gt)
{
	gt->func_tempdata = NULL;
	gt->x_tempdata = NULL;
	gt->rotmatrix = NULL;
	gt->func_shiftxi = NULL;
}

void free_gt(global_t* gt)
{
    if (gt->func_tempdata)	{	free( gt->func_tempdata);	}
    if (gt->x_tempdata)		{	free( gt->x_tempdata );  	}
    if (gt->rotmatrix)		{	 free( gt->rotmatrix);          }
    if (gt->func_shiftxi)	{	free (gt->func_shiftxi);        }
}

/* random_Start(), random_init(), random_exit(), random_Unifo()rm, random_Gauss(), time_tic(), time_tictoc(), time_toc() are adopted 
   from Nikolaus Hansen's source code for CMA-ES
*/

void random_exit(random_t *t)
{
	free( t->rgrand);
}

long random_Start( random_t *t, long unsigned inseed)
{
	long tmp;
	int i;

	t->flgstored = 0;
	t->startseed = inseed;
	if (inseed < 1)
		inseed = 1; 
	t->aktseed = inseed;
#pragma omp parallel for
	for (i = 39; i >= 0; --i)
	{
		tmp = t->aktseed/127773;
		t->aktseed = 16807 * (t->aktseed - tmp * 127773)
			- 2836 * tmp;
		if (t->aktseed < 0) t->aktseed += 2147483647;
		if (i < 32)
			t->rgrand[i] = t->aktseed;
	}
	t->aktrand = t->rgrand[0];
	return inseed;
}

long random_init( random_t *t, long unsigned inseed)
{
	clock_t cloc = clock();

	t->flgstored = 0;
	t->rgrand = (long *)malloc(sizeof(long) * 32);
	if (inseed < 1) {
		while ((long) (cloc - clock()) == 0)
			; 
		inseed = (long)abs(100.0*time(NULL)+clock());
	}
	return random_Start(t, inseed);
}

double random_Uniform( random_t *t)
{
	long tmp;

	tmp = t->aktseed/127773;
	t->aktseed = 16807 * (t->aktseed - tmp * 127773)
		- 2836 * tmp;
	if (t->aktseed < 0) 
		t->aktseed += 2147483647;
	tmp = t->aktrand / 67108865;
	t->aktrand = t->rgrand[tmp];
	t->rgrand[tmp] = t->aktseed;
	return (double)(t->aktrand)/(2.147483647e9);
}

double random_Gauss(random_t *t)
{
	if (USE_ZIGGURAT)
	  return 0;//gsl_ran_gaussian_ziggurat (t);
	double x1, x2, rquad, fac;

	if (t->flgstored)
	{    
		t->flgstored = 0;
		return t->hold;
	}
	do 
	{
		x1 = 2.0 * random_Uniform(t) - 1.0;
		x2 = 2.0 * random_Uniform(t) - 1.0;
		rquad = x1*x1 + x2*x2;
	} while(rquad >= 1 || rquad <= 0);
	fac = sqrt(-2.0*log(rquad)/rquad);
	t->flgstored = 1;
	t->hold = fac * x1;
	return fac * x2;
}

void	time_tic(global_t* t)
{
	t->time_tic_t = time(NULL);	// measure time in seconds
	t->time_tic_c = clock();	// measure time in microseconds up to ~2k seconds
}

double	time_tictoc(global_t* t)
{
	double dt = difftime(t->time_toc_t,t->time_tic_t);
	if (dt < 1000)
		dt = (double)(t->time_toc_c - t->time_tic_c) / CLOCKS_PER_SEC;
	return dt;
}

double	time_toc(global_t* t)
{
	t->time_toc_t = time(NULL);
	t->time_toc_c = clock();
	return time_tictoc(t);
}

void matrix_eye(double* M, int m)
{
    int i, j;
    //if add openmp, slowen the program...
    //#pragma omp parallel for private(j)
	for (i=0; i<m; i++)
		for (j=0; j<m; j++)
			if (i == j)	M[i*m + j] = 1.0;
			else		M[i*m + j] = 0.0;
}


// vector res = matrix a X vector b
void matrix_mult_vector(double* res, double* a, double* b,int m)
{
    int i, j;
	double val = 0.0;
//#pragma omp parallel for private(j)
	for (i=0; i<m; i++)
	{
		val = 0.0;
		for (j=0; j<m; j++)
			val += a[i*m + j] * b[j];
		res[i] = val;
	}
}


// maxtrix res = matrix a X matrix b
void matrix_mult_matrix(double* res, double* a, double* b,int m)
{
    int i, j, k;
	double val;
#pragma omp parallel for private(j), private(k), reduction(+:val)
	for (i=0; i<m; i++)
	    for(j=0; j<m; j++)
		{
			val = 0;
			for (k=0; k<m; k++)
				val += a[i*m + k] * b[k*m + j];
			res[i*m + j] = val;
		}
}


// matrix res = vector a X vector b
void vector_mult_vector(double* res, double* a, double* b,int m)
{
    int i,j;
//#pragma omp parallel for private(j)
	for (i=0; i<m; i++)
		for (j=0; j<m; j++)
			res[i*m + j] = a[i] * b[j];
}


// vector res = vector a X matrix b
void vector_mult_matrix(double* res, double* a, double* b,int m)
{
    int i, j;
	double val;
#pragma omp parallel for private(j), reduction(+:val)
	for (i=0; i<m; i++)
	{
		val = 0;
		for (j=0; j<m; j++)
			val += a[j] * b[j*m + i];
		res[i] = val;
	}
}

double vector_prod(double* a, double* b,int m)
{
    int i;
	double res = 0.0;
//#pragma omp parallel for
	for (i=0; i<m; i++)
		res += a[i] * b[i];
	return res;
}

void generateRotationMatrix(double* B, int N, double* tmp1, random_t* rnd)
{
	double* pB;
	int i, j, k;

//#pragma omp parallel for private(j)
	for (i=0; i<N; i++)
		for (j=0; j<N; j++)
			B[i*N + j] = random_Gauss(rnd);
	
//'without openmp' is faster than 'with openmp'
	for (i=0; i<N; i++)
	{
		for (j=0; j<i; j++)
		{
			double ariarj = 0;
			for (k=0; k<N; k++)
				ariarj = ariarj + B[k*N+i]*B[k*N+j];
			
			for (k=0; k<N; k++)
				B[k*N+i] = B[k*N+i] - ariarj * B[k*N+j];
		}
		double normv = 0;
		for (k=0; k<N; k++)
			normv = normv + B[k*N+i]*B[k*N+i];
			
		normv = sqrt(normv);
		for(k=0; k<N; k++)
			B[k*N+i] = B[k*N+i] / normv;
	}
}

double minv(double a, double b)
{
	if (a < b)	return a;
	else		return b;
}

double maxv(double a, double b)
{
	if (a > b)	return a;
	else		return b;
}

//konohen reduction ayashii
double fsphere(double* x, int N)
{
    int i;
	double Fit = 0;

#pragma omp parallel for reduction(+:Fit)
	for (i=0; i<N; i++)
		Fit += x[i] * x[i];
	return Fit;
}

double felli(double* x, int N)
{
    int i;
	double Fit = 0;
	double alpha = pow(10,6.0);
#pragma omp parallel for reduction(+,Fit)
	for (i=0; i<N; i++)
	    Fit += pow(alpha, (double)i / (double)(N-1) ) * x[i] * x[i];
	return Fit;
}

double felli_fast(double* x, int N, global_t* t)
{
    int i;
	double Fit = 0;
	if (t->func_tempdata == NULL)
	{
	    t->func_tempdata = (double *)malloc(sizeof(double) * N);
		double alpha = pow(10,6.0);
#pragma omp parallel for
		for (i=0; i<N; i++)
		    t->func_tempdata[i] = pow(alpha, (double)i / (double)(N-1) );
	}
#pragma omp parallel for reduction(+:Fit)
	for (i=0; i<N; i++)
		Fit += t->func_tempdata[i] * x[i] * x[i];
	return Fit;
}

double fdiscus(double* x, int N)
{
    int i;
	double Fit = 0;
	Fit = 1e+6 * (x[0] * x[0]);
#pragma omp parallel for reduction(+:Fit)
	for (i=1; i<N; i++)
		Fit += x[i]*x[i];
	return Fit;
}

double fcigar(double* x, int N)
{
    int i;
	double Fit = 0;
#pragma omp parallel for reduction(+:Fit)
	for (i=1; i<N; i++)
		Fit += x[i]*x[i];
	Fit = Fit * 1e+6;
	Fit += x[0] * x[0];
	return Fit;
}

void getRotatedX(double* x, int N, global_t* t)
{
	if (t->x_tempdata == NULL)
	    t->x_tempdata = (double *)malloc(sizeof(double) * N);
	if (t->rotmatrix == NULL)
	{
	    t->rotmatrix = (double *)malloc(sizeof(double) * N*N);
		/*
		FILE* pFile;
		char filename[250];
		sprintf(filename,"matrix%d.txt",N);
		pFile = fopen(filename,"r");
		if (pFile != NULL)
		{
			for (int i = 0; i < N; i++) 
			{
				 for (int j = 0; j < N; j++) 
				 {
					float f;
					fscanf(pFile, "%f", &f);
					t->rotmatrix[i*N + j] = f;
				}
			}
			fclose(pFile);
		}
		else
		{
			for(int i=0; i<N; i++)
				for(int j=0; j<N; j++)
				{
					if (i == j)		t->rotmatrix[i*N + j] = 1;
					else 			t->rotmatrix[i*N + j] = 0;
				}
			printf("CANNOT FIND %s FILE WITH ROTATION MATRIX!!! INDENTITY MATRIX IS USED\n",filename);
		}		
		*/
		generateRotationMatrix(t->rotmatrix, N, t->x_tempdata, &t->ttime);
	}
	matrix_mult_vector(t->x_tempdata, t->rotmatrix, x, N);
}

double frosen(double* x, int N)
{
    int i;
	double Fit = 0;
	double tmp1, tmp2;
	double Fit1 = 0;
	double Fit2 = 0;
	//for (int i=0; i<N-1; i++)
	//	Fit += 100 * pow( x[i]*x[i] - x[i+1], 2.0  ) + pow(x[i] - 1.0, 2.0); // function 'pow' is very slow
	//#pragma omp parallel for reduction(+:Fit1), reduction(+:Fit2)
	//without is faster than with openmp
	for (i=0; i<N-1; i++)
	{
		tmp1 = x[i]*x[i] - x[i+1];
		tmp2 = x[i] - 1.0;
		Fit1 += tmp1*tmp1;
		Fit2 += tmp2*tmp2;
	}
	Fit = 100*Fit1 + Fit2;
	return Fit;
}
int compare(const void * a, const void * b)
{
	if (((sortedvals*)a)->value < ((sortedvals*)b)->value)	return -1;
	if (((sortedvals*)a)->value == ((sortedvals*)b)->value)	return 0;
	if (((sortedvals*)a)->value > ((sortedvals*)b)->value)	return 1;
}

void myqsort(int sz, double* arfitness, int* arindex, sortedvals* arr)
{
    int i;
	for(i=0; i<sz; i++)
	{
		arr[i].value = arfitness[i];
		arr[i].id = i;
	}
	
	qsort( arr, sz, sizeof(sortedvals), compare);
	for(i=0; i<sz; i++)
	{
		arfitness[i] = arr[i].value;
		arindex[i] = arr[i].id;
	}
}


void invAz(int N, double* Av, int iterator_sz, int* iterator, double* v_arr, double* Lj_arr, double K)
{
    int j;
    int p;
	for(j=0; j<iterator_sz; j++)
	{
		int jcur = iterator[j];			
		double* v_j = &v_arr[jcur * N];
		double v_j_mult_Av = 0;
		for(p=0; p<N; p++)					
			v_j_mult_Av += v_j[p] * Av[p];
		v_j_mult_Av = Lj_arr[jcur] * v_j_mult_Av; 
		for(p=0; p<N; p++)
			Av[p] = K * Av[p] - v_j_mult_Av * v_j[p];
	}
}

/*print gene data to file*/
int printGene(FILE *fp, const double *x, int dimension){
    int i;
    for(i=0; i<dimension; ++i){
	fprintf(fp, "%f\t", x[i]);
    }
    return 0;
}/*printGene()*/

void LMCMA(int N, int lambda, int mu, double ccov, double *xmin, double *xmax, int nvectors,
	   int maxsteps, double cc, double val_target, double *sigma, double c_s, double target_f, 
	   int maxevals, int inseed, double* output, int printToFile, int sample_symmetry, MPI_Comm spawn_comm, int num_of_spawn_comm)
{
	// memory allocation
	// m*n
    double* arx = (double *)malloc(sizeof(double) * N * lambda);
    double* v_arr = (double *)malloc(sizeof(double) * N * nvectors);
    double* pc_arr = (double *)malloc(sizeof(double) * N * nvectors);
    // n
    double* pc = (double *)malloc(sizeof(double) * N );
    double* xmean = (double *)malloc(sizeof(double) * N);
    double* xold = (double *)malloc(sizeof(double) * N);
    double* xbestever = (double *)malloc(sizeof(double) * N);
    double* z = (double *)malloc(sizeof(double) * N);
    double* Az = (double *)malloc(sizeof(double) * N);
    double* Av = (double *)malloc(sizeof(double) * N);
    // lambda, mu, nvectors
    double* weights = (double *)malloc(sizeof(double) * mu);
    int* iterator =  (int *)malloc(sizeof(int) * nvectors);		
    double* arfitness = (double *)malloc(sizeof(double) * lambda);
    double* prev_arfitness = (double *)malloc(sizeof(double) * lambda);
    int* arindex = (int *)malloc(sizeof(int) * lambda);
    double* mixed = (double *)malloc(sizeof(double) * lambda * 2);
    int* ranks = (int *)malloc(sizeof(int) * lambda * 2);
    int* ranks_tmp = (int *)malloc(sizeof(int) * lambda * 2);;
    double* Nj_arr = (double *)malloc(sizeof(double) * nvectors);
    double* Lj_arr = (double *)malloc(sizeof(double) * nvectors);
    sortedvals* arr_tmp = (sortedvals *)malloc(sizeof(sortedvals) * 2 * lambda);
    int* t = (int *)malloc(sizeof(int)*nvectors);
    int* vec = (int *)malloc(sizeof(int) * nvectors);
    
    double *scatter_sendvec_w, *scatter_sendvec_d;
    double *scatter_rcvvec_w, *scatter_rcvvec_d;
    double *arFunvals_send, *arFunvals_rcv;
    int num_of_pop_per_spawn = lambda / num_of_spawn_comm;
    double flg_termination = 0;

    my_boundary_transformation_t my_boundaries;
    double *x_temp;
    unsigned int *flg_log;
    flg_log = (unsigned int *)calloc(N, sizeof(unsigned int));

    my_boundary_transformation_init(&my_boundaries, xmin, xmax, flg_log, N);
    x_temp = (double *)calloc(N, sizeof(double));

    scatter_sendvec_w = (double *)malloc(sizeof(double) * (N / 2) * (lambda + num_of_pop_per_spawn));
    scatter_sendvec_d = (double *)malloc(sizeof(double) * (N / 2) * (lambda + num_of_pop_per_spawn));
    arFunvals_send = (double *)malloc(sizeof(double) * num_of_pop_per_spawn);
    arFunvals_rcv = (double *)malloc(sizeof(double) * (lambda + num_of_pop_per_spawn));
    scatter_rcvvec_w = (double *)malloc(sizeof(double) * num_of_pop_per_spawn * N / 2);
    scatter_rcvvec_d = (double *)malloc(sizeof(double) * num_of_pop_per_spawn * N / 2);

    int i, j, k;
    int p;
    global_t gt;
    init_gt(&gt);
	
    // memory initialization
    random_init(&gt.ttime , inseed);

    double sum_weights = 0;
#pragma omp parallel for reduction(+:sum_weights)
    for(i=0; i<mu; i++)
    {
	weights[i] = log((double)(mu+0.5)) - log((double)(1+i));
	sum_weights = sum_weights + weights[i];
    }
    double mueff = 0;
#pragma omp parallel for reduction(+:mueff)
    for(i=0; i<mu; i++)
    {
	weights[i] = weights[i] / sum_weights;
	mueff = mueff + weights[i]*weights[i];
    }
    mueff = 1 / mueff;

    for(i=0; i<N; i++)
	pc[i] = 0;

    double K = 1/sqrt(1 - ccov);
    double M = sqrt(1 - ccov);
	
    for( i=0; i<N; i++)
	xmean[i] = xmin[i] + (xmax[i] - xmin[i])*random_Uniform(&gt.ttime);

    int counteval = 0;
    int iterator_sz = 0;
    double s = 0;	
    int stop = 0;
    int itr = 0;
    int indx = 0;
	
    double BestF;

    FILE* pFile;
    if (printToFile == 1)
    {	
	char filename[250];
	sprintf(filename,"LMCMA%d.txt",N);
	pFile = fopen(filename,"w");
    }
	
      //for test (delete)
    for(i=0;i<((N/2) * (lambda + num_of_pop_per_spawn));++i){
	scatter_sendvec_d[i] = -1;
    }

    while(stop == 0)
    {
	int sign = 1; 
	//this loop is not necessary for the current implementation

	//in main
	for(i=0; i<lambda; ++i){
	    if (sign == 1)
	    { 
		for( k=0; k<N; k++)	// O(n)
		{
		    z[k] = random_Gauss(&gt.ttime);
		    Az[k] = z[k];
		}
			
		for(k=0; k<iterator_sz; k++)	// O(m*n)
		{
		    int jcur = iterator[k];				
		    double* pc_j = &pc_arr[jcur*N];
		    double* v_j = &v_arr[jcur*N];
		    double v_j_mult_z = 0;
		    for(p=0; p<N; p++)	
			v_j_mult_z = v_j_mult_z + v_j[p] * z[p];
		    v_j_mult_z = Nj_arr[jcur] * v_j_mult_z; 				
		    for(p=0; p<N; p++)
			Az[p] = M * Az[p] + v_j_mult_z * pc_j[p];
		}
	    }
	    
	    //generate the gene information (main)
	    //this section was not different from the current implementation largely
	    for(k=0; k<N; k++){	// O(n)
		arx[i*N + k] = xmean[k] + sign*sigma[k]*Az[k];
	    }
	    my_boundary_transformation(&my_boundaries, &arx[i*N], x_temp, 0);
	    for(k=0; k<(int)(N/2); ++k){
		scatter_sendvec_w[i*(N/2)+k+(num_of_pop_per_spawn*(int)(N/2))] = x_temp[k];
		scatter_sendvec_d[i*(N/2)+k+(num_of_pop_per_spawn*(int)(N/2))] = x_temp[k+N/2];
	    }
		    
	    if (sample_symmetry)
		sign = -sign;
	}//end of generate gene information
	
	//printf("scatter information: sendnum = %d\n", num_of_pop_per_spawn * N / 2); 

	//scatter gene information from parent to child, and from child to grandchild	    
	MPI_Scatter(scatter_sendvec_w, num_of_pop_per_spawn * N / 2, MPI_DOUBLE, scatter_rcvvec_w, num_of_pop_per_spawn * N / 2, MPI_DOUBLE, 0, spawn_comm); // in main and make_neuro_spawn
	MPI_Scatter(scatter_sendvec_d, num_of_pop_per_spawn * N / 2, MPI_DOUBLE, scatter_rcvvec_d, num_of_pop_per_spawn * N / 2, MPI_DOUBLE, 0, spawn_comm); // in main and make_neuro_spawn
	//MPI_Scatter(,,,nrn_comm); //in make_neuro_spawn and NEURON
	
	//calculation in NEURON process
	
	//gather fitness information from grandchild to child and from child to parent
	//MPI_Gather(,,,nrn_comm);//in make neuro_spawn and NEURON
	MPI_Gather(arFunvals_send, num_of_pop_per_spawn, MPI_DOUBLE, arFunvals_rcv, num_of_pop_per_spawn, MPI_DOUBLE, 0, spawn_comm);//in main and make_neuro_spawn
		
	for(i=0; i<lambda; ++i){
	    arfitness[i] = arFunvals_rcv[i+num_of_pop_per_spawn];
	    counteval = counteval + 1;
	    if(counteval == 1){
		BestF = arfitness[i];
		for(j=0;j<N;++j){
		    xbestever[j] = arx[i*N+j];
		}
	    }
	    if(BestF > arfitness[i]){
		BestF = arfitness[i];
		for(j=0;j<N;++j){
		    xbestever[j] = arx[i*N+j];
		}
	    }
	}
	
	myqsort(lambda, arfitness, arindex, arr_tmp);

	//save the previous information and reinitialize 'xmean'
	for(i=0; i<N; i++)
	{			
	    xold[i] = xmean[i];
	    xmean[i] = 0;
	}
	//make x_mean
	for(i=0; i<mu; i++)
	{
	    double* cur_x = &arx[arindex[i] * N];
	    for(j=0; j<N; j++)
		xmean[j] += weights[i] * cur_x[j];
	}

	//update pc and Av
	for(i=0; i<N; i++)
	{
	    pc[i] = (1 - cc)*pc[i] + sqrt(cc*(2-cc)*mueff)*(xmean[i] - xold[i])/sigma[i];
	    Av[i] = pc[i];
	}
	//calculate inverse of Av
	invAz(N, Av,iterator_sz, iterator, v_arr, Lj_arr, K);
	
	//printf("itr = %d, nvectors = %d\n", itr, nvectors);
	//??? refer to the original paper
	if (itr < nvectors)
	{
	    t[itr] = itr;	
	}else{
	    int dmin = vec[t[1]] - vec[t[0]];
	    int imin = 1;
	    for(j=1; j<(nvectors-1); j++)
	    {
		int dcur = vec[t[j+1]] - vec[t[j]];
		if (dcur < dmin)
		{
		    dmin = dcur;
		    imin = j + 1;
		}
	    }
	    if (dmin >= maxsteps)
		imin = 0;
	    if (imin != (nvectors-1))
	    {
		int sav = t[imin];
		for(j=imin; j<(nvectors-1); j++)
		    t[j] = t[j+1];
		t[nvectors-1] = sav;
	    }
	}
	
	
	iterator_sz = itr+1;
	if (iterator_sz > nvectors)	iterator_sz = nvectors;
	for(i=0; i<iterator_sz; i++)
	    iterator[i] = t[i];
	int newidx = t[iterator_sz-1];
	vec[newidx] = itr;
			
	for(i=0; i<N; i++)
	{
	    pc_arr[newidx*N + i] = pc[i];
	    v_arr[newidx*N + i] = Av[i];
	}
				
	double nv = 0;
	for(i=0; i<N; i++)
	    nv += Av[i]*Av[i];
	Nj_arr[newidx] = (sqrt(1-ccov)/nv)*(sqrt(1+(ccov/(1-ccov))*nv)-1);
	Lj_arr[newidx] = (1/(sqrt(1-ccov)*nv))*(1-(1/sqrt(1+((ccov)/(1-ccov))*nv)));

	if (itr > 0)
	{
	    for(i=0; i<lambda; i++)
	    {
		mixed[i] = arfitness[i];
		mixed[lambda+i] = prev_arfitness[i];
	    }
	    myqsort(2*lambda, mixed, ranks, arr_tmp);
	    double meanprev = 0;
	    double meancur = 0;
	    for(i=0; i<2*lambda; i++)
		ranks_tmp[i] = ranks[i];
	    for(i=0; i<2*lambda; i++)
		ranks[ranks_tmp[i]] = i;
	    for(i=0; i<lambda; i++)
	    {
		meanprev = meanprev + ranks[i];
		meancur = meancur + ranks[lambda + i];
	    }
	    meanprev = meanprev / lambda;
	    meancur = meancur / lambda;
	    double diffv = (meancur - meanprev)/lambda;
	    double z1 = diffv - val_target;
	    s = (1-c_s)*s + c_s*z1;
	    double d_s = 1;//2.0*(N-1.0)/N;
	    for(i=0;i<N;++i)
		sigma[i] = sigma[i] * exp(s/d_s);
	}

	for(i=0; i<lambda; i++)
	    prev_arfitness[i] = arfitness[i];
				
	if (arfitness[0] < target_f)
	    stop = 1;
	if (counteval >= maxevals)
	    stop = 1;
	itr = itr + 1;
	for(i=0;i<N;++i)	
	    if (sigma[i] < 1e-20)
		stop = 1;
	if ((printToFile == 1) && (pFile))
	    fprintf(pFile,"%d %g\n",counteval,BestF);
	if(counteval >= (lambda * N * 5)){
	    printf("#fbest: %lf\n", BestF);
	    my_boundary_transformation(&my_boundaries, xbestever, x_temp, 0);
	    printGene(stdout, x_temp, N); 
	    printf("\n");
	}
	break;
    }
    flg_termination = 1;
    MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, 0, spawn_comm);
    
    printf("#fbest: %lf\n", BestF);
    my_boundary_transformation(&my_boundaries, xbestever, x_temp, 0);
    printGene(stdout, x_temp, N); 
    printf("\n");
    fflush(stdout);
    my_boundary_transformation_exit(&my_boundaries);

    output[0] = counteval;
    output[1] = BestF;
	
    if (printToFile == 1)
	fclose(pFile);
	
    random_exit(&gt.ttime);
    free(arr_tmp);		free( weights);	free( pc);	free(xmean);	
    free( xold); 			free( z);			free(Az);	free(iterator);	
    free(v_arr);			free(pc_arr);	free(arx);	free( arfitness);
    free(prev_arfitness); free( arindex);	free(Av);	free(t);	
    free(vec);			free(mixed);		free(ranks);	free(ranks_tmp);
    free(Nj_arr);		free(Lj_arr);	
    free_gt(&gt);
    free(scatter_sendvec_w); free(scatter_sendvec_d); free(arFunvals_rcv); free(x_temp);
    free(flg_log); free(xbestever);
}


int loadRangeFile(char *filename, double *xmin_vec, double *xmax_vec, int *N){
    int dimension = 0;
    char buf[TEXT_BUFFER_SIZE];
    double lowerBounds[MAX_NUM_PARAM];
    double upperBounds[MAX_NUM_PARAM];
    unsigned int flg_log[MAX_NUM_PARAM];
    int i;

    FILE *fp;
    if((fp=fopen(filename, "r"))==NULL){
	printf("file open error\n");
	exit(EXIT_FAILURE);
    }
    while( fgets(buf, TEXT_BUFFER_SIZE, fp) != NULL){
	if(strncmp(buf, "#", 1) == 0){ continue;}
	sscanf(buf, "%*s\t%lf\t%lf\t%*lf\t%d\n", &lowerBounds[dimension], &upperBounds[dimension], &flg_log[dimension]);
	//flg_log[dimension] = (unsigned char)atoi((const char*)&flg_log[dimension]);
	dimension++;//make the dimension information here
    }
    *N = dimension;
    for(i=0;i<(*N);++i){
	xmin_vec[i] = lowerBounds[i];
	xmax_vec[i] = upperBounds[i];
    }
    fclose(fp);
    return 0;
}

int main(int argc, char **argv){
    int i;
    int main_myid, main_size;
    char command[] = "./make_neuro_spawn";
    char **spawn_argv;
    int spawn_argv_size=3;
    int num_of_pop_per_child;
    int spawn_numprocs = 6;/*the number of NEURON circuits*/
    MPI_Comm spawn_comm, intercomm, parentcomm;
    int spawn_size, spawn_myid;
    
    int N;                  //      problem dimension
    int lambda;
    int mu;
    double ccov = 1/(10*log((double)N+1.0));//      learning rate for covariance matrix, e.g., 1/(10*log(N+1))
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
    double target_f = 1e-10;        // target fitness function value, e.g., 1e-10
    int maxevals = 1e+8;            // maximum number of function evaluations allowed, e.g., 1e+6
    int inseed = 1;         // initial seed for random number generator, e.g., 1
    int algorithmType = 10; // 0:LMCMA, 10:LMCMA_fixed+features, 1:sepCMA, 2:Cholesky-CMA, 3: baseline CMA-ES
    int printToFile = 0; // 1 or 0

    int sample_symmetry = 0; // use 1 to make the algorithm 2 times faster (in CPU time) and perform better (in terms of number of evaluations)

    double output[2];

    char range_filename[] = "../data/params.txt";

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

    loadRangeFile(range_filename, xmin_tmp, xmax_tmp, &N);
    lambda = 4+floor(3*log((double)N)); //      population size, e.g., 4+floor(3*log(N));
    mu = (int)(lambda / 2); // number of parents, e.g., floor(lambda /2);
    num_of_pop_per_child = lambda / spawn_numprocs;
    nvectors = lambda;
    maxsteps = nvectors;
    cc = (double)(1.0/(double)nvectors);
    printf("N = %d, lambda = %d, mu = %d\n", N, lambda, mu);
    sprintf(spawn_argv[0], "%d", num_of_pop_per_child);
    sprintf(spawn_argv[1], "%d", N);

    xmin_vec = (double *)malloc(sizeof(double) * N);
    xmax_vec = (double *)malloc(sizeof(double) * N);

    for(i=0;i<N;++i){
	xmin_vec[i] = xmin_tmp[i];
	xmax_vec[i] = xmax_tmp[i];
    }
    free(xmin_tmp); xmin_tmp = NULL;
    free(xmax_tmp); xmax_tmp = NULL;

    sigma_vec = (double *)malloc(sizeof(double) * N);
    for(i=0;i<N;++i){
	sigma_vec[i] = 0.5 * (xmax_vec[i] - xmin_vec[i]);
    }

    printf("start comm spawn\n");
    MPI_Comm_spawn(command, spawn_argv, spawn_numprocs, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(intercomm, 0, &spawn_comm);
    MPI_Comm_size(spawn_comm, &spawn_size);
    MPI_Comm_rank(spawn_comm, &spawn_myid);

    printf("start LMCMA\n");
    LMCMA(N, lambda, mu, ccov, xmin_vec, xmax_vec, nvectors, maxsteps, cc, val_target, sigma_vec, c_s, target_f, maxevals, inseed, output, printToFile, sample_symmetry, spawn_comm, spawn_numprocs);
    printf("end of LMCMA\n");

    MPI_Comm_free(&intercomm);
    MPI_Comm_free(&spawn_comm);

    free(xmax_vec); xmax_vec = NULL;
    free(xmin_vec); xmin_vec = NULL;

    MPI_Finalize();
    return 0;
}
