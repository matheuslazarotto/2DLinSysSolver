/* 2D Linear system solver
 * 
 * Reads input data from bash file input.
 * Numerically solve a 2D linear system, 
 *
 * ( \dot{x} )   (A11   A12) ( x )
 * (         ) = (         ) (   )
 * ( \dot{y} )   (A21   A22) ( y )
 *
 * with matrix coefficients [A11, A12, 
 * A21, A22] and return a file 'output.dat' 
 * with the solution.
 * 
 * The plot of the solution is made by a 
 * python script.   
 *
 * Usage: compile with 'make' --> './2dlinsyssolver'
 *  
 * Made by: M. Lazarotto (30/10/2020)
*/

#include <cstring>
#include <chrono>
#include <iomanip>

#include "auxf.hpp"
#include "vector_aux.hpp"
#include "ode_sys_solv.hpp"

void print_prog_bar(int bar_width, double progress, int ASCII_id);

int main(int argc, char **argv) 
{
    /* Start clock */
    auto start = chrono::steady_clock::now();

    /* Parse external argv[] from bash */
    const double a11 = stof(argv[1]);
    const double a12 = stof(argv[2]);
    const double a21 = stof(argv[3]);
    const double a22 = stof(argv[4]);
    const double x_min = stof(argv[5]);
    const double x_max = stof(argv[6]);
    const double y_min = stof(argv[7]);
    const double y_max = stof(argv[8]);
    const double t_run = stof(argv[9]);
    const unsigned int Nx = stoi(argv[10]);
    const unsigned int Ny = stoi(argv[11]);

    /** Integration parameters **/
    double dt0 = 1e-5;                              /* Reference time step */ 
    const double dt_max = 1e-2;                     /* Maximum time step */
    const double dt_print = t_run / (double) 10000; /* Print time frequency */
    const double eps_abs  = 1e-15;                  /* Time step eps absolute control */
    const double eps_rel  = 1e-15;                  /* Time step eps relative control */
    double args[] = {a11, a12, a21, a22};
    const unsigned int dimS = 2;
    
    /** Out file header **/
    FILE *fout = fopen("out_points.dat", "w");      /* Out file: trajectory points */

    fprintf(fout, "; 2D Linear Dynamical System \n");
    fprintf(fout, "; a11 = %f\n", a11);
    fprintf(fout, "; a12 = %f\n", a12);
    fprintf(fout, "; a21 = %f\n", a21);
    fprintf(fout, "; a22 = %f\n", a22);
    fprintf(fout, "; x_min = %f\n", x_min);
    fprintf(fout, "; x_max = %f\n", x_max);
    fprintf(fout, "; y_min = %f\n", y_min);
    fprintf(fout, "; y_max = %f\n", y_max);
    fprintf(fout, "; t_run = %f\n", t_run);

    /** Progress bar settings **/
    std::cout << std::fixed;           /* Set print decimal  */
    std::cout << std::setprecision(1); /* precision.         */
    int Nbar = 20;                     /* Progress bar size  */
    int color_ndx = 31;                /* Progress bar color */

    /** Run **/
    for (int i = 0; i <= (int) Nx; i++)
    {
        for (int j = 0; j <= (int) Ny; j++)
        {
            /* Print progress */
            cout << "\rTrajectory integration (p = ";
            cout << i * j << " / " << Nx * Ny << "): ";
            if (i * j >= (int) (Nx * Ny)) { color_ndx = 32; }
            print_prog_bar(Nbar, (double) (i * j) / (double) (Nx * Ny), color_ndx);
            cout.flush();

            /* Set initial point */
            double xo = x_min + (x_max - x_min) * (double) i / (double) Nx;
            double yo = y_min + (y_max - y_min) * (double) j / (double) Ny;
            double var[] = {xo, yo};
            double t = 0.0;
            double dt = dt0;
            int i_print = 1;

            /** Set GSL system integration **/
            const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rkck;
            gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(type, dimS);
            gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(dimS);
            gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
            gsl_odeiv2_system sys = {system_funct, system_jacob, dimS, &args};

            /* Run dynamics */
            fprintf(fout, "# [t]     [x]     [y]     [dx/dt]    [dy/dt]\n");
            fprintf(fout, "%f    %f    %f    %f    %f\n", t, var[0], var[1], 
                    a11 * var[0] + a12 * var[1], a21 * var[0] + a22 * var[1]);
            while (t <= t_run)
            {
                int status = ode_step_gsl_cntrl(evolve, control, step, sys, 
                                                t, dt, (t+dt_max), var, dimS);
                assert(GSL_SUCCESS == status);

                /* Print trajectory */
                if (t >= i_print * dt_print) 
                {
                    fprintf(fout, "%f    %f    %f    %f    %f\n", t, var[0], var[1], 
                                                        a11 * var[0] + a12 * var[1], 
                                                        a21 * var[0] + a22 * var[1]);

                    i_print += 1;
                }
            }

            gsl_odeiv2_step_free(step);
            gsl_odeiv2_evolve_free(evolve);
            gsl_odeiv2_control_free(control);
        }
    }

    fclose(fout);

    cout << "\nCalculating fixed points and eigenvalues..." << endl;
    
    FILE *fout_eigen = fopen("out_eigen.dat", "w");      /* Out file: eigenstuff */

    /* Fixed points */
    vector<vector<double>> A = {{a11, a12},
                                {a21, a22}};
    vector<double> fixPts(2);
    vector<double> b = {0.0, 0.0};

    lin_sys_solv_LU(A, fixPts, b);

    fprintf(fout_eigen, "# Fixed point: \n");
    fprintf(fout_eigen, "#  fixPt = (%f  %f)\n", fixPts[0], fixPts[1]);
    fprintf(fout_eigen, "#\n");

    /* Eigenvalues & Eigenvectors */
    vector<complex<double>> eigen_vals(2);
    vector<vector<complex<double>>> eigen_vecs(2, vector<complex<double>>(2));

    eigen_stuff_f(A, eigen_vals, eigen_vecs);
    
    fprintf(fout_eigen, "# Eigenvalues: \n");
    for (unsigned int m = 0; m < eigen_vals.size(); m++)
    {
        fprintf(fout_eigen, "#  lamb_%d = %.10f + i * %.10f\n",  
                m, eigen_vals[m].real(), eigen_vals[m].imag());
    }

    fprintf(fout_eigen, "#\n");
    fprintf(fout_eigen, "# Eigenvectors: \n");
    for (unsigned int m = 0; m < eigen_vals.size(); m++) 
    {
        fprintf(fout_eigen, "#  eigen_vec_%d = [", m);
        for (unsigned int n = 0; n < eigen_vals.size(); n++) 
        {
            if (n == (eigen_vals.size() - 1)) 
            {
                fprintf(fout_eigen, "(%lf + i*%lf)]\n", eigen_vecs[n][m].real(), eigen_vecs[n][m].imag());    
            }
            else 
            {
                fprintf(fout_eigen, "(%lf + i*%lf)   ", eigen_vecs[n][m].real(), eigen_vecs[n][m].imag());
            }
        }
    }

    cout << "Done!" << endl;
    
    auto end = chrono::steady_clock::now();
    cout << "\nElapsed time: " << chrono::duration_cast<chrono::minutes>(end - start).count() << "min  "
                               << chrono::duration_cast<chrono::seconds>(end - start).count() % 60 << "sec  "
                               << chrono::duration_cast<chrono::milliseconds>(end - start).count() % 60 << "millisec\n" << endl;

    return 0;
}

void print_prog_bar(int bar_width, double progress, int ASCII_id)
{
    std::cout << "[";
    int pos = (int) (bar_width * progress);
    
    for (int i = 0; i < bar_width; i++)
    {
        if (i <= pos)
        {
            std::cout << "\033[" << ASCII_id << "m";
            std::cout << "=";
            std::cout << "\033[m";
        }
        else
        {
            std::cout << " ";
        }
    }

    std::cout << "] " << int(progress * 100) << " %  \r";
}