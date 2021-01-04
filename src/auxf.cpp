#include "auxf.hpp"

int system_funct(double t, const double s[], double f[], void *args) 
{
    (void)(t); // Avoid unused parameter warning
    double *par = (double *)args;
    double a11 = par[0];
    double a12 = par[1];
    double a21 = par[2];
    double a22 = par[3];
    
    f[0] = a11 * s[0] + a12 * s[1];
    f[1] = a21 * s[0] + a22 * s[1];

    return GSL_SUCCESS; 
}

int system_jacob(double t, const double s[], double *dfdy, double dfdt[], void *args) 
{
    (void)(t); // Avoid unused parameter warning
    double *par = (double *)args;
    double a11 = par[0];
    double a12 = par[1];
    double a21 = par[2];
    double a22 = par[3];

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, a11);
    gsl_matrix_set(m, 0, 1, a12);    
    gsl_matrix_set(m, 1, 0, a21);
    gsl_matrix_set(m, 1, 1, a22);

    dfdt[0] = 0.0; 
    dfdt[1] = 0.0;

    return GSL_SUCCESS;
}