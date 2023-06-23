#include "gsl_debugging.h"

void print_state(int iter, gsl_multiroot_fdfsolver* s)
{
    printf("iter = %3lu x = % .3f % .3f f(x) = % .3e % .3e\n",
        iter,
        gsl_vector_get(s->x, 0),
        gsl_vector_get(s->x, 1),
        gsl_vector_get(s->f, 0),
        gsl_vector_get(s->f, 1));
}
