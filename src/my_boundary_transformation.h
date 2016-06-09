/***************************************
  wrapper of boundary_transformation.c
 ***************************************/

#include "boundary_transformation.h"

#define MAX_NUM_PARAM 150

/*
  CMAES inside values (from -inf to inf, almost from 0 to 10)
  -> temporary values (from 0 to 10)
  -> real values (inside the search range)
 */
typedef struct {
	boundary_transformation_t boundaries;
	double *lower_bounds_real; /* real lower bounds */
	double *upper_bounds_real; /* real upper bounds */
	double *lower_bounds_cmaes; /* real lower bounds */
	double *upper_bounds_cmaes; /* real upper bounds */
        //unsigned char *log_or_not; /* real upper bounds */
        unsigned int *log_or_not; /* real upper bounds */
	unsigned long dimension;
} my_boundary_transformation_t;

void my_boundary_transformation_init(my_boundary_transformation_t *,
									 double const *lower_bounds, double const *upper_bounds, 
									 unsigned int *log_or_not, unsigned long len_of_bounds);

void my_boundary_transformation_exit(my_boundary_transformation_t *);

void my_boundary_transformation(my_boundary_transformation_t *,
				double const *x, double *y, int id); /* new value into y */

/* under construction */
//void my_boundary_transformation_inverse(my_boundary_transformation_t *t,
//		double const *y, double *x, unsigned long len); /* new value into x */
