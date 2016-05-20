#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_boundary_transformation.h"

void my_boundary_transformation_init(my_boundary_transformation_t *t,
									 double const *lower_bounds, double const *upper_bounds, 
									 unsigned char *log_or_not, unsigned long len_of_bounds)
{
	int i;

	t->dimension = len_of_bounds;
	t->lower_bounds_real = (double*)malloc(len_of_bounds * sizeof(double));
	t->upper_bounds_real = (double*)malloc(len_of_bounds * sizeof(double));
	t->lower_bounds_cmaes = (double*)malloc(len_of_bounds * sizeof(double));
	t->upper_bounds_cmaes = (double*)malloc(len_of_bounds * sizeof(double));
	t->log_or_not = (unsigned char*)malloc(len_of_bounds * sizeof(unsigned char));

	for(i = 0; i < len_of_bounds; i++) { 
		t->log_or_not[i] = log_or_not[i];

		if(t->log_or_not[i] != 0) {
			t->lower_bounds_real[i] = log(lower_bounds[i]);
			t->upper_bounds_real[i] = log(upper_bounds[i]);
		} else {
			t->lower_bounds_real[i] = lower_bounds[i];
			t->upper_bounds_real[i] = upper_bounds[i];
		}
		t->lower_bounds_cmaes[i] = 0.0;
		t->upper_bounds_cmaes[i] = 10.0;
	}

	boundary_transformation_init(&t->boundaries, t->lower_bounds_cmaes, t->upper_bounds_cmaes, len_of_bounds);
}

void my_boundary_transformation_exit(my_boundary_transformation_t *t)
{
	boundary_transformation_exit(&t->boundaries);
	free(t->lower_bounds_real);
	free(t->upper_bounds_real);
	free(t->lower_bounds_cmaes);
	free(t->upper_bounds_cmaes);
	free(t->log_or_not);
}

void my_boundary_transformation(my_boundary_transformation_t *t,
								double const *x, double *y)
{
	int i;
	boundary_transformation(&t->boundaries, x, y, t->dimension);
	for(i = 0; i < t->dimension; i++) {
		y[i] = (y[i] - t->lower_bounds_cmaes[i]) 
			* (t->upper_bounds_real[i] - t->lower_bounds_real[i]) 
			/ (t->upper_bounds_cmaes[i] - t->lower_bounds_cmaes[i])
			+ t->lower_bounds_real[i];
		//printf("y[%d] = %lf\n",i,y[i]);
		if(t->log_or_not[i] != 0) { y[i] = exp(y[i]); }
	}
}

/* under construction */
/*
void my_boundary_transformation_inverse(my_boundary_transformation_t *t,
										double const *y, double *x, unsigned long len)
{
}
*/
