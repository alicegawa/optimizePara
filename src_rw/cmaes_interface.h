/* --------------------------------------------------------- */
/* --- File: cmaes_interface.h - Author: Nikolaus Hansen --- */
/* ---------------------- last modified:  IV 2007        --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003, 2007 Nikolaus Hansen. 
     e-mail: hansen AT lri.fr

     License: see file cmaes.c
*/
#include "cmaes.h"

/* --------------------------------------------------------- */
/* ------------------ Interface ---------------------------- */
/* --------------------------------------------------------- */

/* --- initialization, constructors, destructors --- */
double * cmaes_init(cmaes_t *, int dimension , double *xstart, 
					double *stddev, long seed, int lambda, int mu, int maxFunEvals, int maxStopIter,
					const char *input_parameter_filename);
void cmaes_resume_distribution(cmaes_t *evo_ptr, char *filename);
void cmaes_exit(cmaes_t *);

/* --- core functions --- */
double * const * cmaes_SamplePopulation(cmaes_t *);
double *         cmaes_UpdateDistribution(cmaes_t *, 
					  const double *rgFitnessValues);

double *         cmaes_UpdateDistribution_dist(cmaes_t *t, const double *rgFunVal, double *sum_temp, double *sum_cov_temp);

const char *     cmaes_TestForTermination(cmaes_t *);

double *cmaes_SamplePopulation_diag_dist(double *rgD, double sigma, double *x_mean, int dimension, random_t t);
void cmaes_SamplePopulation_diag_dist_update(cmaes_t *t);
random_t cmaes_make_random_t_box(void);

/* --- additional functions --- */
double * const * cmaes_ReSampleSingle( cmaes_t *t, int index);
double const *   cmaes_ReSampleSingle_old(cmaes_t *, double *rgx); 
double *         cmaes_SampleSingleInto( cmaes_t *t, double *rgx);
void             cmaes_UpdateEigensystem(cmaes_t *, int flgforce);

/* --- getter functions --- */
double         cmaes_Get(cmaes_t *, char const *keyword);
const double * cmaes_GetPtr(cmaes_t *, char const *keyword); /* e.g. "xbestever" */
const double * cmaes_GetXPtr(cmaes_t *, unsigned int idx); /* e.g. "xbestever" */
double *       cmaes_GetNew( cmaes_t *t, char const *keyword); /* user is responsible to free */
double *       cmaes_GetInto( cmaes_t *t, char const *keyword, double *mem); /* allocs if mem==NULL, user is responsible to free */

/* --- online control and output --- */
void           cmaes_ReadSignals(cmaes_t *, char const *filename);
void           cmaes_WriteToFile(cmaes_t *, const char *szKeyWord,
                                 const char *output_filename); 
char *         cmaes_SayHello(cmaes_t *);
/* --- misc --- */
double *       cmaes_NewDouble(int n); /* user is responsible to free */
void           cmaes_FATAL(char const *s1, char const *s2, char const *s3, 
			   char const *s4);


