#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#ifdef __cplusplus
  #include "lua.hpp"
#else
  #include "lua.h"
  #include "lualib.h"
  #include "lauxlib.h"
#endif
#include <math.h>

//so that name mangling doesn't mess up function names
#ifdef __cplusplus
extern "C"{
#endif

struct rparams {
    double a;
    double b;
    double c;
    double d;
    double e;
};
  
int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
    printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

int pv_calc_f (const gsl_vector * x, void *params,
              gsl_vector * f)
{
    long double a = ((struct rparams *) params)->a;
    long double b = ((struct rparams *) params)->b;
    long double c = ((struct rparams *) params)->c;
    long double d = ((struct rparams *) params)->d;
    long double e = ((struct rparams *) params)->e;

    const double x0 = gsl_vector_get (x, 0); // i_pv
    const double x1 = gsl_vector_get (x, 1); // i_sh

    const double y0 = a + x0 * b - x1;
    const double y1 = c * exp(x0 * d) + e - x0 - x1;

    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);

    return(GSL_SUCCESS);
}

static int c_pv_calc (lua_State *L) {
    //check and fetch the arguments
    long double arg1 = luaL_checknumber (L, 1);    
    long double arg2 = luaL_checknumber (L, 2);
    long double arg3 = luaL_checknumber (L, 3);
    long double arg4 = luaL_checknumber (L, 4);
    long double arg5 = luaL_checknumber (L, 5);
    long double arg6 = luaL_checknumber (L, 6);
    long double arg7 = luaL_checknumber (L, 7);

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = 2;
    struct rparams p = {arg1, arg2, arg3, arg4, arg5};
    gsl_multiroot_function f = {&pv_calc_f, n, &p};

    long double x_init[2] = {arg6, arg7};
    gsl_vector *x = gsl_vector_alloc (n);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, x);

    print_state (iter, s);

    do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    printf ("status = %s\n", gsl_strerror (status));

    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);

    //push the results
    lua_pushnumber(L, arg2);
    lua_pushnumber(L, arg1);

    //return number of results
    return 2;
}

//library to be registered
static const struct luaL_Reg newlib [] = {
      {"c_pv_calc", c_pv_calc},
      {NULL, NULL}  /* sentinel */
    };

//name of this function is not flexible
int luaopen_newlib (lua_State *L){
    luaL_newlib(L, newlib);
    return 1;
}

#ifdef __cplusplus
}
#endif
