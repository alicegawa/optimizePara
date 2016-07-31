static void Sorted_index(const double *rgFunVal, int *iindex, int n)
{
  int i, j;
  for (i=1, iindex[0]=0; i<n; ++i) {
    for (j=i; j>0; --j) {
      if (rgFunVal[iindex[j-1]] < rgFunVal[i])
        break;
      iindex[j] = iindex[j-1]; /* shift up */
    }
    iindex[j] = i; /* insert i */
  }
}


/* --------------------------------------------------------- */
/* --------------------------------------------------------- */
double *
cmaes_UpdateDistribution_dist( cmaes_t *t, const double *rgFunVal, double *sum_temp)
{
  int i, j, iNk, hsig, N=t->sp.N;
  int flgdiag = ((t->sp.diagonalCov == 1) || (t->sp.diagonalCov >= t->gen)); 
  double sum; 
  double psxps; 
  
  if(t->state == 3)
    FATAL("cmaes_UpdateDistribution(): You need to call \n",
          "SamplePopulation() before update can take place.",0,0);
  if(rgFunVal == NULL) 
    FATAL("cmaes_UpdateDistribution(): ", 
          "Fitness function value array input is missing.",0,0);

  if(t->state == 1)  /* function values are delivered here */
    t->countevals += t->sp.lambda;
  else
    ERRORMESSAGE("cmaes_UpdateDistribution(): unexpected state",0,0,0);
  /* assign function values */
  for (i=0; i < t->sp.lambda; ++i) 
    t->rgrgx[i][N] = t->rgFuncValue[i] = rgFunVal[i];


  /* Generate index */
  Sorted_index(rgFunVal, t->index, t->sp.lambda);
  
  /* Test if function values are identical, escape flat fitness */
  if (t->rgFuncValue[t->index[0]] == 
      t->rgFuncValue[t->index[(int)t->sp.lambda/2]]) {
    t->sigma *= exp(0.2+t->sp.cs/t->sp.damps);
    ERRORMESSAGE("Warning: sigma increased due to equal function values\n",
                 "   Reconsider the formulation of the objective function",0,0);
  }
  
  /* update function value history */
  for(i = (int)*(t->arFuncValueHist-1)-1; i > 0; --i) /* for(i = t->arFuncValueHist[-1]-1; i > 0; --i) */
    t->arFuncValueHist[i] = t->arFuncValueHist[i-1];
  t->arFuncValueHist[0] = rgFunVal[t->index[0]];

  /* update xbestever */
  if (t->rgxbestever[N] > t->rgrgx[t->index[0]][N] || t->gen == 1)
    for (i = 0; i <= N; ++i) {
      t->rgxbestever[i] = t->rgrgx[t->index[0]][i];
      t->rgxbestever[N+1] = t->countevals;
    }

#ifndef _OPENMP
  /* calculate xmean and rgBDz~N(0,C) */
  for (i = 0; i < N; ++i) {
    t->rgxold[i] = t->rgxmean[i]; 
    t->rgxmean[i] = 0.;
    t->rgxmean[i] = sum_temp[i];
    t->rgBDz[i] = sqrt(t->sp.mueff)*(t->rgxmean[i] - t->rgxold[i])/t->sigma; 
  }
#else
  for (i = 0; i < N; ++i) {
	  t->rgxold[i] = t->rgxmean[i];
	  t->rgxmean[i] = sum_temp[i];
	  t->rgBDz[i] = sqrt(t->sp.mueff)*(t->rgxmean[i] - t->rgxold[i])/t->sigma;
  }
#endif /* _OPENMP */

  /* calculate z := D^(-1) * B^(-1) * rgBDz into rgdTmp */
  for (i = 0; i < N; ++i) {
    if (!flgdiag)
      for (j = 0, sum = 0.; j < N; ++j)
        sum += t->B[j][i] * t->rgBDz[j];
    else
      sum = t->rgBDz[i];
    t->rgdTmp[i] = sum / t->rgD[i];
  }
  
  /* cumulation for sigma (ps) using B*z */
  for (i = 0; i < N; ++i) {
    /* if (!flgdiag) */
    /*   for (j = 0, sum = 0.; j < N; ++j) */
    /*     sum += t->B[i][j] * t->rgdTmp[j]; */
    /* else */
      sum = t->rgdTmp[i];
      t->rgps[i] = (1. - t->sp.cs) * t->rgps[i] + 
	  sqrt(t->sp.cs * (2. - t->sp.cs)) * sum;
  }
  
  /* calculate norm(ps)^2 */
  for (i = 0, psxps = 0.; i < N; ++i)
    psxps += t->rgps[i] * t->rgps[i];

  /* cumulation for covariance matrix (pc) using B*D*z~N(0,C) */
  hsig = sqrt(psxps) / sqrt(1. - pow(1.-t->sp.cs, 2*t->gen)) / t->chiN
    < 1.4 + 2./(N+1);
  for (i = 0; i < N; ++i) {
    t->rgpc[i] = (1. - t->sp.ccumcov) * t->rgpc[i] + 
      hsig * sqrt(t->sp.ccumcov * (2. - t->sp.ccumcov)) * t->rgBDz[i];
  }
  
  /* stop initial phase */
  if (t->flgIniphase && 
      t->gen > douMin(1/t->sp.cs, 1+N/t->sp.mucov)) 
    {
      if (psxps / t->sp.damps / (1.-pow((1. - t->sp.cs), t->gen)) 
          < N * 1.05) 
        t->flgIniphase = 0;
    }
    
  /* update of C  */

  Adapt_C2(t, hsig);
  
  /* Adapt_C(t); not used anymore */


  /* update of sigma */
  t->sigma *= exp(((sqrt(psxps)/t->chiN)-1.)*t->sp.cs/t->sp.damps);

  t->state = 3;

  return (t->rgxmean);

} /* cmaes_UpdateDistribution() */

/* --------------------------------------------------------- */
/* --------------------------------------------------------- */
static void
Adapt_C2_dist(cmaes_t *t, int hsig, double *sum_cov_temp)
{
  int i, j, k, N=t->sp.N;
  int flgdiag = ((t->sp.diagonalCov == 1) || (t->sp.diagonalCov >= t->gen)); 
  int count=0;
  
  if (t->sp.ccov != 0. && t->flgIniphase == 0) {
      
      /* definitions for speeding up inner-most loop */
      double ccov1 = douMin(t->sp.ccov * (1./t->sp.mucov) * (flgdiag ? (N+1.5) / 3. : 1.), 1.);
      double ccovmu = douMin(t->sp.ccov * (1-1./t->sp.mucov)* (flgdiag ? (N+1.5) / 3. : 1.), 1.-ccov1); 
      double sigmasquare = t->sigma * t->sigma; 
      
      t->flgEigensysIsUptodate = 0;
      
      /* update covariance matrix */
#ifndef _OPENMP
      for (i = 0; i < N; ++i)
	  for (j = flgdiag ? i : 0; j <= i; ++j) {
	      t->C[i][j] = (1 - ccov1 - ccovmu) * t->C[i][j] 
		  + ccov1
		  * (t->rgpc[i] * t->rgpc[j] 
		     + (1-hsig)*t->sp.ccumcov*(2.-t->sp.ccumcov) * t->C[i][j]);
	      t->C[i][j] = ccovmu * sum_cov_temp[count];
	      ++count;
	  }

#else
      for (i = 0; i < N; ++i) {
	  for (j = flgdiag ? i : 0; j <= i; ++j) {
	      t->C[i][j] = (1 - ccov1 - ccovmu) * t->C[i][j]
		  + ccov1
		  * (t->rgpc[i] * t->rgpc[j]
		     + (1-hsig)*t->sp.ccumcov*(2.-t->sp.ccumcov) * t->C[i][j]);
	      t->C[i][j] = ccovmu * sum_cov_temp[count];
	      ++count;
	  }
      }
#endif /* _OPENMP */
    /* update maximal and minimal diagonal value */
      t->maxdiagC = t->mindiagC = t->C[0][0];
      for (i = 1; i < N; ++i) {
	  if (t->maxdiagC < t->C[i][i])
	      t->maxdiagC = t->C[i][i];
	  else if (t->mindiagC > t->C[i][i])
	      t->mindiagC = t->C[i][i];
      }
  } /* if ccov... */
}
