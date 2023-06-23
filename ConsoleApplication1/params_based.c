#include "params_based.h"

struct rparams {
    double a;
    double b;
    double c;
    double d;
    double e;
};

int pv_calc_f(const gsl_vector* x, void* params,
    gsl_vector* f)
{
    double a = ((struct rparams*)params)->a;
    double b = ((struct rparams*)params)->b;
    double c = ((struct rparams*)params)->c;
    double d = ((struct rparams*)params)->d;
    double e = ((struct rparams*)params)->e;

    const double x0 = gsl_vector_get(x, 0); // i_pv
    const double x1 = gsl_vector_get(x, 1); // i_sh

    const double y0 = a + x0 * b - x1;
    const double y1 = c * exp(x0 * d) + e - x0 - x1;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);

    return(GSL_SUCCESS);
}

int pv_calc_df(const gsl_vector* x, void* params,
    gsl_matrix* J)
{
    const double b = ((struct rparams*)params)->b;
    const double c = ((struct rparams*)params)->c;
    const double d = ((struct rparams*)params)->d;

    const double x0 = gsl_vector_get(x, 0);

    const double df00 = b;
    const double df01 = -1;
    const double df10 = c * d * exp(d * x0) - 1;
    const double df11 = -1;

    gsl_matrix_set(J, 0, 0, df00);
    gsl_matrix_set(J, 0, 1, df01);
    gsl_matrix_set(J, 1, 0, df10);
    gsl_matrix_set(J, 1, 1, df11);

    return GSL_SUCCESS;
}

int pv_calc_fdf(const gsl_vector* x, void* params,
    gsl_vector* f, gsl_matrix* J)
{
    pv_calc_f(x, params, f);
    pv_calc_df(x, params, J);

    return GSL_SUCCESS;
}

void c_pv_calc(double arg1, double arg2, double arg3,
    double arg4, double arg5, double arg6, double arg7,
    double* ret_1, double* ret_2) {

    const gsl_multiroot_fdfsolver_type* T;
    gsl_multiroot_fdfsolver* s;

    int status;

    const size_t n = 2;
    struct rparams p = { arg1, arg2, arg3, arg4, arg5 };
    gsl_multiroot_function_fdf f = { &pv_calc_f,
                                    &pv_calc_df,
                                    &pv_calc_fdf,
                                    n, &p };

    double x_init[2] = { arg6, arg7 };
    gsl_vector* x = gsl_vector_alloc(n);

    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fdfsolver_gnewton;
    s = gsl_multiroot_fdfsolver_alloc(T, n);
    gsl_multiroot_fdfsolver_set(s, &f, x);
    int iter = 0;

    //print_state(iter, s);

    do
    {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(s);

        //print_state(iter, s);

        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);

    //printf("status = %s\n", gsl_strerror(status));

    *ret_1 = gsl_vector_get(s->x, 0);
    *ret_2 = gsl_vector_get(s->x, 1);

    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);
}

void solve_params_based() {
    // Constant values
    const double q = 1.6e-19; // elementary charge
    const double K = 1.4e-23; // Boltzmann constant
    const double T_n = 298; // Kelvin

    // Datasheet values (inputs)
    double n = 1.3; // Diode factor
    double N_s = 60; // Amount of PV cells in series
    double N_p = 1; // Amount of PV cells in parallel
    double i_sc = 8.42; // short circuit current(A)
    double k_i = 0.0032; // short circuit thermal gain
    double v_oc = 36.6; // open circuit voltage of PV module(V)
    double R_s = 0.221; // series resistance(Ohm)
    double R_sh = 415.405; // shunt resistance(Ohm)
    double E_go = 1.1; // bandgap semiconductor(eV)
    double T = 25; // degrees C
    double G = 1000; // solar irradiance

    // Calculated values
    T = T + 273.15;
    // photocurrent
    double i_ph = (i_sc + k_i * (T - 298)) * G / 1000.0;
    /// reverse saturation current
    double i_rs = i_sc / (exp(q * v_oc / (n * N_s * K * T)) - 1);
    // saturation current
    double i_0 = i_rs * pow((T / T_n), 3) * exp(q * E_go * (1 / T_n - 1 / T) / (n * K));

    // Array declaration
    double i[100] = { 0 };
    double p[100] = { 0 };
    double v[100] = { 0 };

    v_oc = v_oc / N_s;
    double b = R_s / R_sh;
    double d = q * R_s / (n * K * N_s * T);
    double e = i_0 + i_ph;
    double i_pv, i_sh;
    // Generate values for a range of the equation
    int index = 0;
    double v_pv = 0;
    for (v_pv; v_pv < 52; v_pv += v_oc) {
        double a = v_pv / R_sh;
        double c = -i_0 * exp(q * v_pv / (n * K * N_s * T));
        c_pv_calc(a, b, c, d, e, 0, 0, &i_pv, &i_sh);
        if (i_pv >= 0) {
            i[index] = i_pv;
            p[index] = i_pv * v_pv;
            v[index] = v_pv;
            index++;
        }
    }

    for (int iter = 0; iter < index; iter++) {
        printf("%f ", i[iter]);
    }

    for (int iter = 0; iter < index; iter++) {
        printf("%f ", p[iter]);
    }

    for (int iter = 0; iter < index; iter++) {
        printf("%f ", v[iter]);
    }
}
