#ifndef _ODE_SYS_H_
#define _ODE_SYS_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;

int ode_step_gsl_cntrl(gsl_odeiv2_evolve *evolve, gsl_odeiv2_control *control, 
                       gsl_odeiv2_step *step, gsl_odeiv2_system &sys, double &t, 
                       double &dt, double t_max, double *s, unsigned int dimS);

int ode_step_gsl_fixed(gsl_odeiv2_driver *driver, double &t, double dt, double *s,
                       unsigned int dimS, int nSteps);

#endif