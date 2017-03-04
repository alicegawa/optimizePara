#include "mplmcma.h"

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

int printGene2(FILE *fp, const double *x, int dimension){
    int i;
    for(i=0; i<dimension; ++i){
	fprintf(fp, "%f\n", x[i]);
    }
    return 0;
}/*printGene2()*/

void nullset2arrays(mplmcma_t *t){
  t->arx=NULL;
  t->v_arr=NULL;
  t->pc_arr=NULL;
  t->pc=NULL;
  t->xmean=NULL;
  t->xold=NULL;
  t->xbestever=NULL;
  t->z=NULL;
  t->Az=NULL;
  t->Av=NULL;

  t->weights=NULL;
  t->iterator=NULL;
  t->arfitness=NULL;
  t->prev_arfitness=NULL;
  t->arindex=NULL;
  t->mixed=NULL;
  t->ranks=NULL;
  t->ranks_tmp=NULL;
  t->Nj_arr=NULL;
  t->Lj_arr=NULL;
  t->arr_tmp=NULL;
  t->ti=NULL;
  t->vec=NULL;
   
  /*for communicate with the child processt*/
  t->scatter_sendvec_w, t->scatter_sendvec_d=NULL;
  t->scatter_sendvec=NULL;
  t->scatter_rcvvec_w=NULL;
  t->scatter_rcvvec_d=NULL;
  t->scatter_rcvvec=NULL;
  t->arFunvals_send=NULL;
  t->arFunvals_rcv=NULL;

  t->x_temp=NULL;
  t->flg_log=NULL;
}
int init_lmcma_arrays(mplmcma_t *t){
  nullset2arrays(t);
  // m*n
  t->arx = (double *)malloc(sizeof(double) * t->N * t->lambda);
  t->v_arr = (double *)malloc(sizeof(double) * t->N * t->nvectors);
  t->pc_arr = (double *)malloc(sizeof(double) * t->N * t->nvectors);
  // n
  t->pc = (double *)malloc(sizeof(double) * t->N );
  t->xmean = (double *)malloc(sizeof(double) * t->N);
  t->xold = (double *)malloc(sizeof(double) * t->N);
  t->xbestever = (double *)malloc(sizeof(double) * t->N);
  t->z = (double *)malloc(sizeof(double) * t->N);
  t->Az = (double *)malloc(sizeof(double) * t->N);
  t->Av = (double *)malloc(sizeof(double) * t->N);
  // t->lambda, mu, t->nvectors
  t->weights = (double *)malloc(sizeof(double) * t->mu);
  t->iterator =  (int *)malloc(sizeof(int) * t->nvectors);		
  t->arfitness = (double *)malloc(sizeof(double) * t->lambda);
  t->prev_arfitness = (double *)malloc(sizeof(double) * t->lambda);
  t->arindex = (int *)malloc(sizeof(int) * t->lambda);
  t->mixed = (double *)malloc(sizeof(double) * t->lambda * 2);
  t->ranks = (int *)malloc(sizeof(int) * t->lambda * 2);
  t->ranks_tmp = (int *)malloc(sizeof(int) * t->lambda * 2);;
  t->Nj_arr = (double *)malloc(sizeof(double) * t->nvectors);
  t->Lj_arr = (double *)malloc(sizeof(double) * t->nvectors);
  t->arr_tmp = (sortedvals *)malloc(sizeof(sortedvals) * 2 * t->lambda);
  t->ti = (int *)malloc(sizeof(int)*t->nvectors);
  t->vec = (int *)malloc(sizeof(int) * t->nvectors);
    
  t->num_of_pop_per_spawn = t->lambda / t->spawn_numprocs;//num_of_spawn_comm;
  t->flg_termination = 0;

  t->flg_log = (unsigned int *)calloc(t->N, sizeof(unsigned int));
  t->x_temp = (double *)calloc(t->N, sizeof(double));

  t->scatter_sendvec_w = (double *)malloc(sizeof(double) * (t->N / 2) * (t->lambda + t->num_of_pop_per_spawn));
  t->scatter_sendvec_d = (double *)malloc(sizeof(double) * (t->N / 2) * (t->lambda + t->num_of_pop_per_spawn));
  t->scatter_sendvec = (double *)malloc(sizeof(double) * t->N * (t->lambda + t->num_of_pop_per_spawn));
  t->arFunvals_send = (double *)malloc(sizeof(double) * t->num_of_pop_per_spawn);
  t->arFunvals_rcv = (double *)malloc(sizeof(double) * (t->lambda + t->num_of_pop_per_spawn));
  t->scatter_rcvvec_w = (double *)malloc(sizeof(double) * t->num_of_pop_per_spawn * t->N / 2);
  t->scatter_rcvvec_d = (double *)malloc(sizeof(double) * t->num_of_pop_per_spawn * t->N / 2);
  t->scatter_rcvvec = (double *)malloc(sizeof(double) * t->num_of_pop_per_spawn * t->N);
  if(t->arx==NULL || t->v_arr==NULL || t->pc_arr==NULL || t->pc==NULL || t->xmean==NULL || t->xold==NULL || t->xbestever==NULL || t->z==NULL || t->Az==NULL || t->Av==NULL || t->weights==NULL || t->iterator==NULL || t->arfitness==NULL || t->prev_arfitness==NULL || t->arindex==NULL || t->mixed==NULL || t->ranks==NULL || t->ranks_tmp==NULL || t->Nj_arr==NULL || t->Lj_arr==NULL || t->arr_tmp==NULL || t->ti==NULL || t->vec==NULL || t->flg_log==NULL || t->x_temp==NULL || t->scatter_sendvec_w==NULL || t->scatter_sendvec_d==NULL || t->scatter_sendvec==NULL || t->arFunvals_send==NULL || t->arFunvals_rcv==NULL || t->scatter_rcvvec_w==NULL || t->scatter_rcvvec_d==NULL || t->scatter_rcvvec==NULL){
    printf("memory allocation error occurs in init_lmcma_arrays\n");
    return 0;
  }else{
    return 1;
  }
}

void free_lmcma_arrays(mplmcma_t *t){
  free(t->arx);	
  free(t->v_arr);
  free(t->pc_arr);
  
  free(t->pc);
  free(t->xmean); 
  free(t->xold);
  free(t->xbestever);
  free(t->z);
  free(t->Az);
  free(t->Av);
  free(t->weights);
  free(t->iterator);
  free(t->arfitness);
  free(t->prev_arfitness); 
  free(t->arindex);
  free(t->mixed);
  free(t->ranks);
  free(t->ranks_tmp);
  free(t->Nj_arr);
  free(t->Lj_arr);
  free(t->arr_tmp);	
  free(t->ti);
  free(t->vec);

  free(t->x_temp);
  free(t->flg_log);
  
  free(t->scatter_sendvec_w); 
  free(t->scatter_sendvec_d); 
  free(t->scatter_sendvec);
  free(t->arFunvals_send);
  free(t->arFunvals_rcv); 
  free(t->scatter_rcvvec_w); 
  free(t->scatter_rcvvec_d); 
  free(t->scatter_rcvvec); 

  nullset2arrays(t);
}
//void LMCMA(int N, int lambda, int mu, double ccov, double *xmin, double *xmax, int nvectors,
//	   int maxsteps, double cc, double val_target, double *sigma, double c_s, double target_f, 
//	   int maxevals, int inseed, double* output, int printToFile, int sample_symmetry, MPI_Comm spawn_comm, int num_of_spawn_comm,// mplmcma_t *mplmcma)
void LMCMA(mplmcma_t *t)
{
  // memory allocation
  int res;
  res = init_lmcma_arrays(t);
  if(!res){
    return ;
  }
    my_boundary_transformation_t my_boundaries;
    my_boundary_transformation_init(&my_boundaries, t->xmin, t->xmax, t->flg_log, t->N);
    
    int i, j, k;
    int p;
    global_t gt;
    init_gt(&gt);
	
    // memory initialization
    random_init(&gt.ttime , t->inseed);

    double sum_weights = 0;
#pragma omp parallel for reduction(+:sum_weights)
    for(i=0; i<t->mu; i++)
    {
	t->weights[i] = log((double)(t->mu+0.5)) - log((double)(1+i));
	sum_weights = sum_weights + t->weights[i];
    }
    double mueff = 0;
#pragma omp parallel for reduction(+:mueff)
    for(i=0; i<t->mu; i++)
    {
	t->weights[i] = t->weights[i] / sum_weights;
	mueff = mueff + t->weights[i]*t->weights[i];
    }
    mueff = 1 / mueff;

    for(i=0; i<t->N; i++)
	t->pc[i] = 0;

    double K = 1/sqrt(1 - t->ccov);
    double M = sqrt(1 - t->ccov);
	
    for( i=0; i<t->N; i++)
	t->xmean[i] = t->xmin[i] + (t->xmax[i] - t->xmin[i])*random_Uniform(&gt.ttime);

    int counteval = 0;
    int iterator_sz = 0;
    double s = 0;	
    int stop = 0;
    int itr = 0;
    int indx = 0;
	
    double BestF=100000;

    FILE* pFile;
    if (t->printToFile == 1)
    {	
	char filename[250];
	sprintf(filename,"LMCMA%d.txt",t->N);
	pFile = fopen(filename,"w");
    }
	
      //for test (delete)
    for(i=0;i<((t->N/2) * (t->lambda + t->num_of_pop_per_spawn));++i){
	t->scatter_sendvec_d[i] = -1;
    }

    int loop_count = 0;
    while(stop == 0)
    {
	int sign = 1; 
	//this loop is not necessary for the current implementation

	//in main
	for(i=0; i<t->lambda; ++i){
	    if (sign == 1)
	    { 
		for( k=0; k<t->N; k++)	// O(n)
		{
		    t->z[k] = random_Gauss(&gt.ttime);
		    t->Az[k] = t->z[k];
		}
			
		for(k=0; k<iterator_sz; k++)	// O(m*n)
		{
		    int jcur = t->iterator[k];				
		    double* pc_j = &t->pc_arr[jcur*t->N];
		    double* v_j = &t->v_arr[jcur*t->N];
		    double v_j_mult_z = 0;
		    for(p=0; p<t->N; p++){
			v_j_mult_z = v_j_mult_z + v_j[p] * t->z[p];
		    }
		    v_j_mult_z = t->Nj_arr[jcur] * v_j_mult_z;
		    //printf("\n\n\nt->Nj_arr[%d] = %lf\n\n\n\n", jcur, t->Nj_arr[jcur]);
		    for(p=0; p<t->N; p++)
			t->Az[p] = M * t->Az[p] + v_j_mult_z * pc_j[p];
		}
	    }
	    
	    //generate the gene information (main)
	    //this section was not different from the current implementation largely
	    for(k=0; k<t->N; k++){	// O(n)
		t->arx[i*t->N + k] = t->xmean[k] + sign*t->sigma[k]*t->Az[k];
		//		printf("t->arx[%d] = %lf\n", i*t->N+k, t->arx[i*t->N+k]);
	    }
	    my_boundary_transformation(&my_boundaries, &t->arx[i*t->N], t->x_temp, 0);
	    /* for(k=0; k<(int)(t->N/2); ++k){ */
	    /* 	t->scatter_sendvec_w[i*(t->N/2)+k+(t->num_of_pop_per_spawn*(int)(t->N/2))] = t->x_temp[k]; */
	    /* 	t->scatter_sendvec_d[i*(t->N/2)+k+(t->num_of_pop_per_spawn*(int)(t->N/2))] = t->x_temp[k+t->N/2]; */
	    /* } */
	    for(k=0; k<t->N; ++k){
		t->scatter_sendvec[i*t->N+k+(t->num_of_pop_per_spawn * t->N)] = t->x_temp[k];
	    }
		    
	    if (t->sample_symmetry)
		sign = -sign;
	}//end of generate gene information
	
	//printf("t->scatter information: sendnum = %d\n", t->num_of_pop_per_spawn * t->N / 2); 

	//scatter gene information from parent to child, and from child to grandchild	    
//	MPI_Scatter(t->scatter_sendvec_w, t->num_of_pop_per_spawn * t->N / 2, MPI_DOUBLE, t->scatter_rcvvec_w, t->num_of_pop_per_spawn * t->N / 2, MPI_DOUBLE, 0, t->spawn_comm); // in main and make_neuro_spawn
//	MPI_Scatter(t->scatter_sendvec_d, t->num_of_pop_per_spawn * t->N / 2, MPI_DOUBLE, t->scatter_rcvvec_d, t->num_of_pop_per_spawn * t->N / 2, MPI_DOUBLE, 0, t->spawn_comm); // in main and make_neuro_spawn
	MPI_Scatter(t->scatter_sendvec, t->num_of_pop_per_spawn * t->N, MPI_DOUBLE, t->scatter_rcvvec, t->num_of_pop_per_spawn * t->N, MPI_DOUBLE, 0, t->spawn_comm);
	//MPI_Scatter(,,,nrn_comm); //in make_neuro_spawn and NEURON
	
	//calculation in NEURON process
	
	//gather fitness information from grandchild to child and from child to parent
	//MPI_Gather(,,,nrn_comm);//in make neuro_spawn and NEURON
	MPI_Gather(t->arFunvals_send, t->num_of_pop_per_spawn, MPI_DOUBLE, t->arFunvals_rcv, t->num_of_pop_per_spawn, MPI_DOUBLE, 0, t->spawn_comm);//in main and make_neuro_spawn
		
	for(i=0; i<t->lambda; ++i){
	    t->arfitness[i] = t->arFunvals_rcv[i+t->num_of_pop_per_spawn];
	    counteval = counteval + 1;
	    if(counteval == 1){
		BestF = t->arfitness[i];
		for(j=0;j<t->N;++j){
		    t->xbestever[j] = t->arx[i*t->N+j];
		}
	    }
	    if(BestF > t->arfitness[i]){
		BestF = t->arfitness[i];
		for(j=0;j<t->N;++j){
		    t->xbestever[j] = t->arx[i*t->N+j];
		}
	    }
	}
	
	myqsort(t->lambda, t->arfitness, t->arindex, t->arr_tmp);

	//save the previous information and reinitialize 'xmean'
	for(i=0; i<t->N; i++)
	{			
	    t->xold[i] = t->xmean[i];
	    t->xmean[i] = 0;
	}
	//make x_mean
	for(i=0; i<t->mu; i++)
	{
	    double* cur_x = &t->arx[t->arindex[i] * t->N];
	    for(j=0; j<t->N; j++)
		t->xmean[j] += t->weights[i] * cur_x[j];
	}

	//update pc and Av
	for(i=0; i<t->N; i++)
	{
	    t->pc[i] = (1 - t->cc)*t->pc[i] + sqrt(t->cc*(2-t->cc)*mueff)*(t->xmean[i] - t->xold[i])/t->sigma[i];
	    t->Av[i] = t->pc[i];
	}
	//calculate inverse of Av
	invAz(t->N, t->Av,iterator_sz, t->iterator, t->v_arr, t->Lj_arr, K);
	
	//printf("itr = %d, t->nvectors = %d\n", itr, t->nvectors);
	//??? refer to the original paper
	if (itr < t->nvectors)
	{
	    t->ti[itr] = itr;	
	}else{
	    int dmin = t->vec[t->ti[1]] - t->vec[t->ti[0]];
	    int imin = 1;
	    for(j=1; j<(t->nvectors-1); j++)
	    {
		int dcur = t->vec[t->ti[j+1]] - t->vec[t->ti[j]];
		if (dcur < dmin)
		{
		    dmin = dcur;
		    imin = j + 1;
		}
	    }
	    if (dmin >= t->maxsteps)
		imin = 0;
	    if (imin != (t->nvectors-1))
	    {
		int sav = t->ti[imin];
		for(j=imin; j<(t->nvectors-1); j++)
		    t->ti[j] = t->ti[j+1];
		t->ti[t->nvectors-1] = sav;
	    }
	}
	
	
	iterator_sz = itr+1;
	if (iterator_sz > t->nvectors)	iterator_sz = t->nvectors;
	for(i=0; i<iterator_sz; i++)
	    t->iterator[i] = t->ti[i];
	int newidx = t->ti[iterator_sz-1];
	t->vec[newidx] = itr;
			
	for(i=0; i<t->N; i++)
	{
	    t->pc_arr[newidx*t->N + i] = t->pc[i];
	    t->v_arr[newidx*t->N + i] = t->Av[i];
	}
				
	double nv = 0;
	for(i=0; i<t->N; i++)
	    nv += t->Av[i]*t->Av[i];
	//printf("1-t->ccov = %lf, nv = %lf, 1+(t->ccov/(1-t->ccov))*nv = %lf\n", 1-t->ccov, nv, 1+(t->ccov/(1-t->ccov))*nv); 
	t->Nj_arr[newidx] = (sqrt(1-t->ccov)/nv)*(sqrt(1+(t->ccov/(1-t->ccov))*nv)-1);
	t->Lj_arr[newidx] = (1/(sqrt(1-t->ccov)*nv))*(1-(1/sqrt(1+((t->ccov)/(1-t->ccov))*nv)));

	if (itr > 0)
	{
	    for(i=0; i<t->lambda; i++)
	    {
		t->mixed[i] = t->arfitness[i];
		t->mixed[t->lambda+i] = t->prev_arfitness[i];
	    }
	    myqsort(2*t->lambda, t->mixed, t->ranks, t->arr_tmp);
	    double meanprev = 0;
	    double meancur = 0;
	    for(i=0; i<2*t->lambda; i++)
		t->ranks_tmp[i] = t->ranks[i];
	    for(i=0; i<2*t->lambda; i++)
		t->ranks[t->ranks_tmp[i]] = i;
	    for(i=0; i<t->lambda; i++)
	    {
		meanprev = meanprev + t->ranks[i];
		meancur = meancur + t->ranks[t->lambda + i];
	    }
	    meanprev = meanprev / t->lambda;
	    meancur = meancur / t->lambda;
	    double diffv = (meancur - meanprev)/t->lambda;
	    double z1 = diffv - t->val_target;
	    s = (1-t->c_s)*s + t->c_s*z1;
	    double d_s = 1;//2.0*(t->N-1.0)/t->N;
	    for(i=0;i<t->N;++i)
		t->sigma[i] = t->sigma[i] * exp(s/d_s);
	}

	for(i=0; i<t->lambda; i++)
	    t->prev_arfitness[i] = t->arfitness[i];
				
	if (t->arfitness[0] < t->target_f){
	  stop = 1;
	  printf("arfitness is small adequately\n");
	}
	if (counteval >= t->maxevals){
	    stop = 1;
	    printf("evaluation is max times\n");
	}
	itr = itr + 1;
	for(i=0;i<t->N;++i){	
	  if (t->sigma[i] < 1e-20){
	    stop = 1;
	    printf("sigma is small adequately\n");
	  }
	  if(t->sigma[i] > 10000000){
	    t->sigma[i] = 80000;
	  }
	  if(loop_count%100==0){
	    t->sigma[i] = 10.0;
	  }
	}
	if ((t->printToFile == 1) && (pFile))
	    fprintf(pFile,"%d %g\n",counteval,BestF);
	/* if(counteval >= (t->lambda * t->N * 5)){ */
	/*     printf("#fbest: %lf\n", BestF); */
	/*     my_boundary_transformation(&my_boundaries, t->xbestever, t->x_temp, 0); */
	/*     printGene(stdout, t->x_temp, t->N);  */
	/*     printf("\n"); */
	/* } */
	if(stop==1){
	    t->flg_termination = 1;
	}
	loop_count++;
	printf("%d times loop finish. fbest = %lf\n", loop_count, BestF);
	//t->flg_termination=1;
	MPI_Bcast(&t->flg_termination, 1, MPI_DOUBLE, 0, t->spawn_comm);
	//break;
    }
    //t->flg_termination = 1;

    
    printf("#fbest: %lf\n", BestF);
    my_boundary_transformation(&my_boundaries, t->xbestever, t->x_temp, 0);
    printGene(stdout, t->x_temp, t->N);
    FILE *fp_result;
    if((fp_result=fopen("result.txt", "w"))==NULL){
      printf("file open error occurs @fp_result\n");
    }
    printGene2(fp_result, t->x_temp, t->N);
    fclose(fp_result);
    printf("\n");
    fflush(stdout);
    my_boundary_transformation_exit(&my_boundaries);

    t->output[0] = counteval;
    t->output[1] = BestF;
	
    if (t->printToFile == 1)
	fclose(pFile);
	
    random_exit(&gt.ttime);
    free_gt(&gt);
    free_lmcma_arrays(t);
}


int loadRangeFile(char *filename, double *xmin_vec, double *xmax_vec, int *N){
    int dimension = 0;
    char buf[TEXT_BUFFER_SIZE];
    double lowerBounds[MAX_NUM_PARAM];
    double upperBounds[MAX_NUM_PARAM];
    unsigned int flg_log[MAX_NUM_PARAM];
    int i;

    FILE *fp=NULL;
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
	printf("lowerBounds[%d] = %lf, upperBounds[%d] = %lf, flg_log[%d] = %d\n", i, lowerBounds[i], i, upperBounds[i], i, flg_log[i]);
    }
    fclose(fp);
    return 0;
}

int init_hyperparam_arrays(mplmcma_t *t){
  t->xmin = (double *)malloc(sizeof(double) * t->N);
  t->xmax = (double *)malloc(sizeof(double) * t->N);
  t->sigma = (double *)malloc(sizeof(double) * t->N);
}

void free_hyperparam_arrays(mplmcma_t *t){
  free(t->xmin); t->xmin = NULL;
  free(t->xmax); t->xmax = NULL;
  free(t->sigma); t->sigma = NULL;
}
