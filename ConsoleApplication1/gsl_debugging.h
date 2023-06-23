#ifndef GSL_DEBUGGING_H
#define GSL_DEBUGGING_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <math.h>

void print_state(int iter, gsl_multiroot_fdfsolver* s);

#endif
