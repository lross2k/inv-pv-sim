#include "single_diode.h"

void solve_single_diode() {
    // Constant values
    const double q = 1.6e-19; // elementary charge
    const double k = 1.4e-23; // Boltzmann constant
    const double T_r = 25 + 273.15; // Kelvin

    // Datasheet values (inputs)
    double A = 1.2; // Ideality factor
    double n_p = 60; // Amount of PV cells in PV module
    double i_sc = 8.49; // short circuit current (A)
    double alpha = 6; // mA / C
    double T = 25 + 273.15; // degrees C
    double S = 1000; // solar irradiance
    double v_oc = 37.5; // open circuit voltage of PV module(V)

    // Calculated values
    v_oc = v_oc / n_p; // open circuit voltage of PV Cell(V);
    double i_rr = exp(log(i_sc) - v_oc * q / (A * k * T));
    double E_g = 1.16 - 7.02e-4 * (T * T) / (T - 1108);
    double i_rs = i_rr * pow((T / T_r), 3) * exp(q * E_g / (k * A) * (1 / T_r - 1 / T));
    double n_s = n_p;
    double i_ph = (i_sc + alpha * (T - T_r)) * S / 1e3; //Photovoltaic current source

    // Array declaration
    double i[100] = { 0 };
    double p[100] = { 0 };
    double v[100] = { 0 };

    // Generate values for a range of the equation
    int index = 0;
    double v_p = 0;
    printf("%f < %f\n", v_p, v_oc);
    for (v_p; v_p < 52; v_p += v_oc) {
        double i_calc = i_ph - i_rs * (exp(q * v_p / (A * k * T * n_s)) - 1);
        if (i_calc >= 0) {
            i[index] = i_calc;
            p[index] = i_calc * v_p;
            v[index] = v_p;
            index = index + 1;
        } else {
            break;
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
