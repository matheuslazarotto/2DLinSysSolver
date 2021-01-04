#include "ode_sys_solv.hpp"

int ode_step_gsl_cntrl(gsl_odeiv2_evolve *evolve, gsl_odeiv2_control *control, 
                       gsl_odeiv2_step *step, gsl_odeiv2_system &sys, 
                       double &t, double &dt, double t_max, double *s, 
                       unsigned int dimS)
{
    /** performs ode step from [t] to [t_max] with precision control and 
     *  return success status **/    
    return gsl_odeiv2_evolve_apply(evolve, control, step, 
                                   &sys, &t, t_max, &dt, s);
}

int ode_step_gsl_fixed(gsl_odeiv2_driver *driver, double &t, double dt, double *s,
                       unsigned int dimS, int nSteps) 
{
    /** performs [nSteps] ODE steps with fixed time step and return success 
     *  status **/    
    return gsl_odeiv2_driver_apply_fixed_step(driver, &t, dt, nSteps, s);

}