#ifndef _AUXF_H_
#define _AUXF_H_

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

using namespace std;

/** Dynamical system functions **/
int system_funct(double t, const double s[], double f[], void *args);
int system_jacob(double t, const double s[], double *dfdy, double dfdt[], void *args);

#endif