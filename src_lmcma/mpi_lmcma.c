#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

void LMCMA(int N, int lambda, int mu, double ccov, double *xmin, double *xmax, int nvectors,
	int maxsteps, double cc, double val_target, double sigma, double c_s, double target_f, 
	   int maxevals, enum FunctionId FuncId, int inseed, double* output, int printToFile, int sample_symmetry, MPI_Comm spawn_comm, int num_of_spawn_comm)
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

    scatter_sendvec_w = (double *)malloc(sizeof(double) * (N / 2) * (lambda + 1));
    scatter_sendvec_d = (double *)malloc(sizeof(double) * (N / 2) * (lambda + 1));
    arFunvals_rcv = (double *)malloc(sizeof(double) * lambda * (num_of_spawn_comm + 1));
  
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
	sprintf(filename,"LMCMA%dfunc%d.txt",N,(int)FuncId);
	pFile = fopen(filename,"w");
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
	    for(k=0; k<N; k++)	// O(n)
		arx[i*N + k] = xmean[k] + sign*sigma*Az[k];
	    for(k=0; k<(int)(N/2); ++k){
		scatter_rcvvec_w[i*(N/2)+k] = arx[i*N + k];
		scatter_rcvvec_d[i*(N/2)+k] = arx[i*N + k + N / 2];
	    }
		    
	    if (sample_symmetry)
		sign = -sign;
	}//end of generate gene information

	//scatter gene information from parent to child, and from child to grandchild	    
	MPI_Scatter(scatter_sendvec_w, num_of_spawn_comm * N / 2, MPI_DOUBLE, scatter_rcvvec_w, num_of_spawn_comm * N / 2, MPI_DOUBLE, 0, spawn_comm); // in main and make_neuro_spawn
	MPI_Scatter(scatter_sendvec_w, num_of_spawn_comm * N / 2, MPI_DOUBLE, scatter_rcvvec_d, num_of_spawn_comm * N / 2, MPI_DOUBLE, 0, spawn_comm); // in main and make_neuro_spawn
	//MPI_Scatter(,,,nrn_comm); //in make_neuro_spawn and NEURON

	//calculation in NEURON process
	
	//gather fitness information from grandchild to child and from child to parent
	//MPI_Gather(,,,nrn_comm);//in make neuro_spawn and NEURON
	MPI_Gather(arFunvals_send, num_of_pop_per_spawn, MPI_DOUBLE, arFunvals_rcv, num_of_pop_per_spawn, MPI_DOUBLE, 0, spawn_comm);//in main and make_neuro_spawn
	
	for(i=0; i<lambda; ++i){
	    arfitness[i] = arFunvals_rcv[i+num_of_pop_per_spawn];
	    counteval = counteval + 1
	    if(counteval == 1){
		BestF = arfitness[i];
	    }
	    if(BestF > arfitness[i]){
		BestF = arfitness[i];
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
	    pc[i] = (1 - cc)*pc[i] + sqrt(cc*(2-cc)*mueff)*(xmean[i] - xold[i])/sigma;
	    Av[i] = pc[i];
	}
	//calculate inverse of Av
	invAz(N, Av,iterator_sz, iterator, v_arr, Lj_arr, K);


	//??? refer to the original paper
	if (itr < nvectors)
	{
	    t[itr] = itr;			
	}
	else
	{
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
	    sigma = sigma * exp(s/d_s);
	}

	for(i=0; i<lambda; i++)
	    prev_arfitness[i] = arfitness[i];
				
	if (arfitness[0] < target_f)
	    stop = 1;
	if (counteval >= maxevals)
	    stop = 1;
	itr = itr + 1;
		
	if (sigma < 1e-20)
	    stop = 1;
	if ((printToFile == 1) && (pFile))
	    fprintf(pFile,"%d %g\n",counteval,BestF);
    }
	
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
    free(scatter_sendvec_w); free(scatter_sendvec_d); free(arFunvals_rcv);
}
