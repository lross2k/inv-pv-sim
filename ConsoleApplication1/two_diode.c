#include "two_diode.h"

struct td_rparams {
    double a, b, c, d, e, f, g, h;
};

double td_pv_calc_f(double x, void* params)
{
    const double a = ((struct td_rparams*)params)->a;
    const double b = ((struct td_rparams*)params)->b;
    const double c = ((struct td_rparams*)params)->c;
    const double d = ((struct td_rparams*)params)->d;
    const double e = ((struct td_rparams*)params)->e;
    const double pf = ((struct td_rparams*)params)->f;
    const double g = ((struct td_rparams*)params)->g;
    const double h = ((struct td_rparams*)params)->h;

    return(a - b * (exp((pf + x * g) / (c * e)) - 1) - 
           b * (exp((pf + x * g) / (d * e)) - 1) - (pf + x * g) / h - x);
}

double td_pv_calc_df(double x, void* params)
{
    const double a = ((struct td_rparams*)params)->a;
    const double b = ((struct td_rparams*)params)->b;
    const double c = ((struct td_rparams*)params)->c;
    const double d = ((struct td_rparams*)params)->d;
    const double e = ((struct td_rparams*)params)->e;
    const double f = ((struct td_rparams*)params)->f;
    const double g = ((struct td_rparams*)params)->g;
    const double h = ((struct td_rparams*)params)->h;

    return(- b * exp(f / (c * e)) * exp(x * g / (c * e)) * g / (c * e) -
           b * exp(f / (d * e)) * exp(x * g / (d * e)) * g / (d * e) + g / h - 1);
}

void td_pv_calc_fdf(double x, void* params,
                   double *y, double *dy)
{
    const double a = ((struct td_rparams*)params)->a;
    const double b = ((struct td_rparams*)params)->b;
    const double c = ((struct td_rparams*)params)->c;
    const double d = ((struct td_rparams*)params)->d;
    const double e = ((struct td_rparams*)params)->e;
    const double pf = ((struct td_rparams*)params)->f;
    const double g = ((struct td_rparams*)params)->g;
    const double h = ((struct td_rparams*)params)->h;

    *y = a - b * (exp((pf + x * g) / (c * e)) - 1) -
         b * (exp((pf + x * g) / (d * e)) - 1) - (pf + x * g) / h - x;
    *dy = -b * exp(pf / (c * e)) * exp(x * g / (c * e)) * g / (c * e) -
        b * exp(pf / (d * e)) * exp(x * g / (d * e)) * g / (d * e) + g / h - 1;
}

void solve_two_diode_current(double I_pv, double I_01, double a_1,
                             double a_2, double V_T, double V, double R_s, 
                             double R_p, double I, double* ret) {

    int status;
    int iter = 0, max_iter = 10000;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver* s;
    double x0, x = 0;
    gsl_function_fdf FDF;
    struct td_rparams params = { I_pv, I_01, a_1, a_2, V_T, V, R_s, R_p };

    FDF.f = &td_pv_calc_f;
    FDF.df = &td_pv_calc_df;
    FDF.fdf = &td_pv_calc_fdf;
    FDF.params = &params;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_root_fdfsolver_set(s, &FDF, x);

    //printf("using %s method\n",
    //    gsl_root_fsolver_name(s));

    //printf("%5s [%9s, %9s] %9s\n",
    //    "iter", "lower", "upper", "root");

    do
    {
        iter++;
        status = gsl_root_fdfsolver_iterate(s);
        x0 = x;
        x = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(x, x0, 0, 1e-3);

        //if (status == GSL_SUCCESS)
            //printf("Converged:\n");

        //printf("%5d %10.7f\n",
        //    iter, x);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    *ret = x;

    gsl_root_fdfsolver_free(s);
}

void solve_two_diode() {
    // Constant values
    const double q = 1.6e-19; // elementary charge
    const double K = 1.4e-23; // Boltzmann constant

    // Datasheet values (inputs)
    double G_n = 1000;
    double T_n = 298; // Kelvin
    double a_1 = 1; //ideality factor diode(assumed)
    double a_2 = 2; //ideality factor diode 2 (assumed)
    double N_s = 36; //Amount of PV cells in series
    double i_scn = 3.8; // short circuit current(A)
    double k_i = 0.003; // short circuit thermal gain(A / K)
    double k_v = -80e-3;  // (V / K)
    double v_oc = 21.1; //open circuit voltage of PV module(V)
    double i_mp = 3.5; //maximum power point current(A)
    double v_mp = 17.1; //maximum power point voltage(V)
    double p_max = 60; //maximum power(W)

    //input values
    double T = 25; //degrees C
    double G = 1000; //solar irradiance

    // Calculated values
    double v_oc_mp = v_oc;
    v_oc = v_oc / N_s;
    T = T + 273.15; //Kelvin
    double V_T = N_s * K * T / q;
    double R_s = 0;
    double I_pvn = 1 / i_mp; // Starting value equals to 1 / i_mp
    double DeltaT = T - T_n;
    double I_0 = (i_mp + k_i * DeltaT) / exp((v_oc_mp + k_v * DeltaT) / (a_1 * V_T) - 1);
    double I_01 = I_0 * (i_mp + k_i * DeltaT) / exp(v_mp + k_v * DeltaT / ((a_1 + a_2) * V_T) - 1);
    double P_maxe = i_mp * v_mp;
    double I_pv = (I_pvn + k_i * DeltaT) * (G / G_n);
    double R_p = v_mp + i_mp * R_s / (I_pv - I_0 * (exp(v_mp + i_mp * R_s / V_T) + exp(v_mp + i_mp * R_s / V_T) - 2));
    double I_sc = (i_scn + k_i * DeltaT) * (G / G_n);
    //Solve the system with newton rhapson

    // Array declaration
    double i[100] = { 0 };
    double p[100] = { 0 };
    double v[100] = { 0 };

    // Generate values for a range of the equation
    int index = 0;
    double I = 0;
    for (double V = 0; V < 52; V += v_oc) {
        solve_two_diode_current(I_pv, I_01, a_1, a_2, V_T, V, R_s, R_p, I, &I);
        I = I_pv - I_01 * (exp((V + I * R_s) / (a_1 * V_T)) - 1) - I_01 * (exp((V + I * R_s) / (a_2 * V_T)) - 1) - (V + I * R_s) / R_p;
        //if (I >= 0) {
            i[index] = I;
            p[index] = I * V;
            v[index] = V;
            index++;
        //}
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

    FILE* f = NULL;
    fopen_s(&f, "two_diode.txt", "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(f, "i = [");
    for (int iter = 0; iter < index - 1; iter++) {
        fprintf(f, "%.5f, ", i[iter]);
    }
    fprintf(f, "%.5f]\n", i[index - 1]);

    fprintf(f, "p = [");
    for (int iter = 0; iter < index - 1; iter++) {
        fprintf(f, "%.5f, ", p[iter]);
    }
    fprintf(f, "%.5f]\n", p[index - 1]);

    fprintf(f, "v = [");
    for (int iter = 0; iter < index - 1; iter++) {
        fprintf(f, "%.5f, ", v[iter]);
    }
    fprintf(f, "%.5f]\n", v[index - 1]);

    fclose(f);
}
